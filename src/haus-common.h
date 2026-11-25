#ifndef SMILE_HAUS_COMMON_H_
#define SMILE_HAUS_COMMON_H_

#include "geom-common.h"
#include <string>

// Input conversion, part splitting, and 2D polygon primitives shared by
// src/haus-r2.cpp (in planar coordinates) and src/haus-s2.cpp (in gnomonic
// tangent-plane coordinates).

namespace smile_haus {

using smile_geom::Pt2;

enum HausMethod {
  kDiscrete = 0,
  kBoundary = 1,
  kPolygon = 2
};

inline HausMethod parse_haus_method(const std::string& method, double tol) {
  HausMethod code = kDiscrete;
  if (method == "discrete") {
    code = kDiscrete;
  } else if (method == "boundary") {
    code = kBoundary;
  } else if (method == "polygon") {
    code = kPolygon;
  } else {
    Rcpp::stop(
      "Unknown method: must be \"discrete\", \"boundary\", or \"polygon\"."
    );
  }
  if (code != kDiscrete && (!std::isfinite(tol) || tol <= 0.0)) {
    Rcpp::stop("tol must be a single finite positive number.");
  }
  return code;
}

// Converts an R object to a numeric matrix with at least 2 columns, or stops.
inline Rcpp::NumericMatrix as_unit_matrix(SEXP x) {
  Rcpp::NumericMatrix m(x);
  if (m.nrow() > 0 && m.ncol() < 2) {
    Rcpp::stop("Each unit must be a matrix with at least 2 columns.");
  }
  return m;
}

// Returns a zero-copy Eigen view of `m`.
inline Eigen::Map<const Eigen::MatrixXd> as_eigen(const Rcpp::NumericMatrix& m) {
  return Eigen::Map<const Eigen::MatrixXd>(m.begin(), m.nrow(), m.ncol());
}

template <typename UnitT, typename MakeFn>
std::vector<UnitT> make_unit_vec(const Rcpp::List& polys, MakeFn make_fn) {
  std::vector<UnitT> out;
  out.reserve(polys.size());
  for (R_xlen_t i = 0; i < polys.size(); i++) {
    const Rcpp::NumericMatrix m = as_unit_matrix(polys[i]);
    out.push_back(make_fn(as_eigen(m)));
  }
  return out;
}

inline bool has_valid_coords(const Eigen::Ref<const Eigen::MatrixXd>& poly) {
  const int nr = static_cast<int>(poly.rows());
  if (nr == 0) return false;
  for (int i = 0; i < nr; i++) {
    if (!std::isfinite(poly(i, 0)) || !std::isfinite(poly(i, 1))) return false;
  }
  return true;
}

struct RawPart {
  int begin;
  int end;
  bool is_poly;
  bool close_ring;
};

// Splits a unit matrix into contiguous parts (rings, lines, or points).
inline std::vector<RawPart> extract_raw_parts(
    const Eigen::Ref<const Eigen::MatrixXd>& poly) {
  std::vector<RawPart> parts;
  const int nr = static_cast<int>(poly.rows());
  if (nr == 0) return parts;
  const bool has_part = poly.cols() > 2;
  const bool has_poly_col = poly.cols() > 3;

  int p_start = 0;
  while (p_start < nr) {
    int p_end = p_start + 1;
    if (!has_part) {
      p_end = nr;
    } else {
      const double pid = poly(p_start, 2);
      while (p_end < nr && poly(p_end, 2) == pid) p_end++;
    }
    const int plen = p_end - p_start;
    const bool is_poly =
      has_poly_col ? (poly(p_start, 3) != 0.0) : (plen >= 3);
    const bool close_ring =
      !has_part && plen > 1 &&
      (poly(0, 0) != poly(nr - 1, 0) || poly(0, 1) != poly(nr - 1, 1));
    parts.push_back({p_start, p_end, is_poly, close_ring});
    p_start = p_end;
  }
  return parts;
}

struct PolyPart2D {
  int e_begin;
  int e_end;
  bool is_poly;
  double xmin, ymin, xmax, ymax;
};

// Canonical unit ordering so continuous Hausdorff distances are bit-for-bit
// symmetric between (x, y) and (y, x).
template <typename KeyT>
inline bool should_swap_units(int nv_x, int ne_x, const KeyT& key_x,
                              int nv_y, int ne_y, const KeyT& key_y) {
  if (nv_y != nv_x) return nv_y < nv_x;
  if (ne_y != ne_x) return ne_y < ne_x;
  return key_y < key_x;
}

// Exact 2D segment vs. axis-aligned box intersection (Liang-Barsky), shifted
// relative to the box center to preserve precision on large coordinates.
inline bool seg_intersects_box(double ax, double ay, double bx, double by,
                               double x0, double y0, double x1, double y1) {
  if (std::max(ax, bx) < x0 || std::min(ax, bx) > x1 ||
      std::max(ay, by) < y0 || std::min(ay, by) > y1) {
    return false;
  }
  const double cx = 0.5 * (x0 + x1);
  const double cy = 0.5 * (y0 + y1);
  const double hx = 0.5 * (x1 - x0);
  const double hy = 0.5 * (y1 - y0);
  const double px = ax - cx;
  const double py = ay - cy;
  const double dx = bx - ax;
  const double dy = by - ay;
  double t0 = 0.0;
  double t1 = 1.0;

  auto clip = [&](double p, double q) -> bool {
    if (p == 0.0) return q >= 0.0;
    const double r = q / p;
    if (p < 0.0) {
      if (r > t1) return false;
      if (r > t0) t0 = r;
    } else {
      if (r < t0) return false;
      if (r < t1) t1 = r;
    }
    return true;
  };

  return clip(-dx, px + hx) && clip(dx, hx - px) &&
         clip(-dy, py + hy) && clip(dy, hy - py) && (t0 <= t1);
}

// Even-odd ray casting over all polygon rings of a unit. Shifts coordinates by
// (px, py) before multiplying to avoid cancellation.
inline bool point_in_poly_2d(double px, double py,
                             const std::vector<Pt2>& pts,
                             const std::vector<int>& edge_a,
                             const std::vector<int>& edge_b,
                             const std::vector<PolyPart2D>& parts,
                             double xmin, double ymin,
                             double xmax, double ymax,
                             double pad = 0.0) {
  if (px < xmin - pad || px > xmax + pad ||
      py < ymin - pad || py > ymax + pad) {
    return false;
  }
  bool inside = false;
  const int np = static_cast<int>(parts.size());
  for (int p = 0; p < np; p++) {
    const PolyPart2D& part = parts[p];
    if (!part.is_poly) continue;
    if (px < part.xmin - pad || px > part.xmax + pad ||
        py < part.ymin - pad || py > part.ymax + pad) {
      continue;
    }
    for (int k = part.e_begin; k < part.e_end; k++) {
      const int ib = edge_b[k];
      if (ib < 0) continue;
      const int ia = edge_a[k];
      const double y1 = pts[ia].y - py;
      const double y2 = pts[ib].y - py;
      if ((y1 > 0.0) != (y2 > 0.0)) {
        const double x1 = pts[ia].x - px;
        const double x2 = pts[ib].x - px;
        const double x_int = x1 - y1 * (x2 - x1) / (y2 - y1);
        if (x_int > 0.0) inside = !inside;
      }
    }
  }
  return inside;
}

}  // namespace smile_haus

#endif  // SMILE_HAUS_COMMON_H_
