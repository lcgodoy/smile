## Makefile for generating R packages.
## Based on Rob J Hyndman's Makefile:
##   https://github.com/robjhyndman/forecast/blob/master/Makefile

PKG_NAME=$(shell grep -i ^package DESCRIPTION | cut -d : -d \  -f 2)

.PHONY: render_md build_site build install clean docs check check_cran

default: build

check:
	Rscript -e "devtools::check()"

check_cran:
	R CMD BUILD .
	R CMD check --as-cran $(PKG_NAME)_*.tar.gz
	rm -f .$(PKG_NAME)_*.tar.gz

render_md:
	Rscript -e "knitr::knit('README.Rmd')"

build_site:
	Rscript -e "pkgdown::build_site()"

build:
	Rscript -e "devtools::build()"

install:
	R CMD INSTALL .

clean:
	rm -f ./$(PKG_NAME)_*.tar.gz
	rm -rf ./$(PKG_NAME).Rcheck
	rm -rf man/*.Rd
	rm -rf NAMESPACE
	rm -f *~ R/*~ src/*~
	rm -f R/*~ src/*~
	rm -f src/*~

docs:
	Rscript -e "devtools::document()"

# revdep:
# 	Rscript -e "revdepcheck::revdep_check(num_workers=3)"

release:
	Rscript -e "devtools::submit_cran(args = c('--compact-vignettes=both'))"
