# 4. Method

### General setup

Let us set up the notations first. Suppose a there exists a partition of
a region $D \in \mathcal{R}^{2}$ (e.g., a city). This partition is
denoted by $A_{i}$, $i = 1,\ldots,n$. Moreover, there exists another
partition of the same city, denoted $B_{j}$, where $j = 1,\ldots,m$.
These partitions can be seen as two different administrative divisions
within the same city. It is common for different government agencies to
release data for different divisions of a same city, country, or state.

### Model-based approach

Consequently,

$$E\left\lbrack Y\left( A_{i} \right) \right\rbrack = \frac{1}{|A_{i}|}\int_{A_{i}}E\left\lbrack Y(\mathbf{s}) \right\rbrack\, d\mathbf{s} = \frac{1}{|A_{i}|}\int_{A_{i}}\mu\, d\mathbf{s} = \mu$$

and

$$\begin{aligned}
{{Cov}\left\lbrack Y\left( A_{i} \right),Y\left( A_{j} \right) \right\rbrack} & {= \frac{1}{|A_{i}||A_{j}|}\int_{A_{i} \times A_{j}}{Cov}\left\lbrack Y(\mathbf{s},Y(\mathbf{s}\prime) \right\rbrack\, d\mathbf{s}\, d{\mathbf{s}\mathbf{\prime}}} \\
 & {= \frac{1}{|A_{i}||A_{j}|}\int_{A_{i} \times A_{j}}C( \parallel \mathbf{s} - \mathbf{s}\prime \parallel ;{\mathbf{θ}})\, d\mathbf{s}\, d{\mathbf{s}\mathbf{\prime}},}
\end{aligned}$$

where $\parallel \mathbf{s} - \mathbf{s}\prime \parallel$ is the
Euclidean distance between the coordinates $\mathbf{s}$ and
$\mathbf{s}\prime$, and
$C( \parallel \mathbf{s} - \mathbf{s}\prime \parallel ;{\mathbf{θ}})$ is
an isotropic covariance function depending on the parameter
$\mathbf{θ}$.

Assume we observe a random variable $Y( \cdot )$ at each region $A_{i}$
and we are interested in predict/estimate this variable in each of the
regions $B_{j}$. Now suppose the random variable $Y( \cdot )$ varies
continuously over $D$ and is defined as follows
$$Y(\mathbf{s}) = \mu + S(\mathbf{s}) + \varepsilon(\mathbf{s}),\,\mathbf{s} \in D \subset \mathcal{R}^{2}.$$

where
$$S( \cdot ) \sim {GP}\left( 0,\sigma^{2}\rho( \cdot ;\,\phi,\kappa) \right)\;{\mspace{6mu}\text{and}\mspace{6mu}}\;\varepsilon( \cdot )\overset{i.i.d.}{\sim}N\left( 0,\sigma^{2}\rho( \cdot ;\,\phi,\kappa) \right),$$
where $S$ and $\varepsilon$ are independent. For now, let’s make the
unrealistic assumption that all those parameters are known. Then, our
assumption is that the observed data is as follows

$$\begin{aligned}
{Y\left( A_{i} \right)} & {= \frac{1}{|A_{i}|}\int_{A_{i}}Y(\mathbf{s})\, d\mathbf{s}} \\
 & {= \frac{1}{|A_{i}|}\int_{A_{i}}\left\lbrack \mu + S(\mathbf{s}) + \varepsilon(\mathbf{s}) \right\rbrack\, d\mathbf{s}} \\
 & {= \mu + \frac{1}{|A_{i}|}\int_{A_{i}}S(\mathbf{s})d\mathbf{s} + \frac{1}{|A_{i}|}\int_{A_{i}}\varepsilon(\mathbf{s})d\mathbf{s},}
\end{aligned}$$

where $| \cdot |$ returns the area of a polygon. Furthermore, it can be
shown that (using Fubini’s Theorem and some algebraic manipulation)
$${Cov}\left( Y\left( A_{i} \right),Y\left( A_{j} \right) \right) = \frac{\sigma^{2}}{|A_{i}||A_{j}|}\int_{A_{i} \times A_{j}}\rho( \parallel \mathbf{s} - \mathbf{s}\prime \parallel ;\,\phi,\kappa)\, d\mathbf{s}\, d\mathbf{s}\prime + \mathbf{I}(i = j)\frac{\tau}{|A_{i}|},$$
where $\rho( \cdot ;\,\phi,\kappa)$ is a positive definite correlation
function. Now, let $R_{\kappa}(\phi)$ be a correlation matrix such that
$$R_{\kappa}(\phi)_{ij} = \frac{1}{|A_{i}||A_{j}|}\int_{A_{i} \times A_{j}}\rho( \parallel \mathbf{s} - \mathbf{s}\prime \parallel ;\,\phi,\kappa)\, d\mathbf{s}\, d\mathbf{s}\prime,$$
thus,
$$Y\left( A_{1},\cdots,A_{n} \right) \sim N\left( \mu\mathbf{1}_{n},\sigma^{2}R_{\kappa}(\phi) + \tau{diag}\left( |A_{1}|^{- 1},\ldots,|A_{1}|^{- 1} \right) \right).$$
Then, if we assume
$\left( Y^{\top}\left( A_{1},\cdots,A_{n} \right),Y^{\top}\left( B_{1},\cdots,A_{m} \right)^{\top} \right)$
to be jointly normal, we use can the conditional mean of
$Y^{\top}\left( B_{1},\cdots,A_{m} \right)^{\top}$ given
$Y^{\top}\left( A_{1},\cdots,A_{n} \right)$ to estimate the observed
random variable in the partition $B_{1},\ldots,B_{m}$.

------------------------------------------------------------------------

Now, suppose the parameters
${\mathbf{θ}} = \left( \mu,\sigma^{2},\phi,\tau \right)$ are unknown.
The Likelihood of $Y\left( A_{1},\ldots,A_{n} \right)$ can still be
computed.

In particular, if we use the parametrization $\alpha = \tau/\sigma^{2}$,
we have closed form for the Maximum Likelihood estimators both for $\mu$
and $\sigma^{2}$. Thus, we can optimize the profile likelihood for
$\phi$ and $\alpha$ numerically. Then, we resort on conditional Normal
properties again to compute the predictions in a new different set of
regions.

### Areal Interpolation (AI)

Areal interpolation is a nonparametric approach that interpolates
$Y\left( A_{i} \right)$’s to construct $Y\left( B_{j} \right)$’s. Define
an $m \times n$ matrix $\mathbf{W} = \{ w_{ij}\}$, where $w_{ij}$ is the
weight associated with the polygon $A_{i}$ in constructing
$Y\left( B_{j} \right)$. The weights are
$w_{ij} = |A_{i} \cap B_{j}|/|B_{j}|$(Goodchild and Lam 1980; Gotway and
Young 2002). The interpolation for
$\widehat{Y}\left( B_{1},\ldots,B_{m} \right)$ is constructed as
$$\widehat{Y}\left( B_{1},\ldots,B_{m} \right) = \mathbf{W}Y\left( A_{1},\ldots,A_{n} \right).$$
The expectation and variance of the predictor are, respectively,
$$E\left\lbrack \widehat{Y}\left( B_{1},\ldots,B_{m} \right) \right\rbrack = \mathbf{W}E\left\lbrack Y\left( A_{1},\ldots,A_{n} \right) \right\rbrack$$
and
$$\text{Var}\left\lbrack \widehat{Y}\left( B_{1},\ldots,B_{m} \right) \right\rbrack = \mathbf{W}\text{Var}\left\lbrack Y\left( A_{1},\ldots,A_{n} \right) \right\rbrack\mathbf{W}^{\top}.$$
In practice, the covariance matrix
$\text{Var}\left\lbrack Y\left( A_{1},\ldots,A_{n} \right) \right\rbrack$
is unknown and, consequently needs to be estimated.

The variance each predictor
$\text{Var}\left\lbrack \widehat{Y}\left( B_{i} \right) \right\rbrack$
is needed as an uncertainty measure. It relies on both the variances of
$Y\left( A_{j} \right)$’s and their covariances: $$\begin{array}{r}
{\text{Var}\left\lbrack \widehat{Y}\left( B_{i} \right) \right\rbrack = \sum\limits_{i = 1}^{n}w_{ij}^{2}\text{Var}\left\lbrack Y\left( A_{i} \right) \right\rbrack + 2\sum\limits_{l \neq i}w_{ij}w_{il}\text{Cov}\left\lbrack Y\left( A_{i} \right),Y\left( A_{l} \right) \right\rbrack.}
\end{array}$$ The variances are often observed in survey data, but the
covariances are not. For practical purpose, we propose an approximation
for
$\text{Cov}\left\lbrack Y\left( A_{i} \right),Y\left( A_{l} \right) \right\rbrack$
based on Moran’s I, a global spatial autocorrelation. Specifically, let
$\rho_{I}$ be the Moran’s I calculated with a weight matrix constructed
with first-degree neighbors. That is, $\rho_{I}$ is the average of the
pairwise correlation for all neighboring pairs. For regions $A_{i}$ and
$A_{l}$, if they are neighbors of each other, our approximation is
$$\begin{array}{r}
{\text{Cov}\left\lbrack Y\left( A_{i} \right),Y\left( A_{l} \right) \right\rbrack = \rho_{I}\sqrt{\text{Var}\left\lbrack Y\left( A_{i} \right) \right\rbrack\text{Var}\left\lbrack Y\left( A_{l} \right) \right\rbrack}.}
\end{array}$$ The covariance between non-neighboring $A_{i}$ and $A_{l}$
is discarded. The final uncertainty approximation for
$\text{Var}\left\lbrack \widehat{Y}\left( B_{i} \right) \right\rbrack$
will be an underestimate. Alternatively, we can derive, at least, an
upper bound for the variance of the estimates by using a simple
application from the Cauchy–Schwartz inequality, in which case,
$\rho_{I}$ is replaced with~1.

## Reference

Goodchild, Michael F, and Nina Siu-Ngan Lam. 1980. “Areal Interpolation:
A Variant of the Traditional Spatial Problem.” *Geo-Processing* 1:
279–312.

Gotway, Carol A, and Linda J Young. 2002. “Combining Incompatible
Spatial Data.” *Journal of the American Statistical Association* 97
(458): 632–48.
