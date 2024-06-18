include config.mk

### * --- Development

## ----PKG DEVELOPMENT-----

### ** todo

## todo : Search 'TODO' (in '*.[R|Rmd|Rmd.original|stan|stanfunctions]')
.PHONY: todo
todo:
	find . -name "*.R" -o -name "*.Rmd" -o -name "*.Rmd.original" -o -name "*.stan" -o -name "*.stanfunctions" | xargs grep -i TODO --color=always

### ** recompile

## recompile : recompile Stan code and install pkg (no "data", no "document")
.PHONY: recompile
recompile:
	@$(MAKE) clean-compiled
	@$(MAKE) onlyinstall

### ** test

## test : run pkg tests
.PHONY: test
test:
	@printf "\n"
	@printf "$(GREEN)=== Running package tests ===$(NC)\n"
	@printf "\n"
	@Rscript -e "r <- testthat::ProgressReporter[[\"new\"]](show_praise = FALSE); devtools::test(reporter = r)"
	@cd tests/testthat; rm -f Rplots.pdf

### ** test-examples

## examples : run pkg examples for testing (also "document")
.PHONY: examples
examples: document
	@printf "\n"
	@printf "$(GREEN)=== Running package examples ===$(NC)\n"
	@printf "\n"
	@Rscript -e "devtools::run_examples()"
	@rm -f Rplots.pdf

### ** test-cs

## test-cs : test case studies
.PHONY: test-cs
test-cs:
	@printf "\n"
	@printf "$(GREEN)*** Testing case studies ***$(NC)\n"
	@printf "\n"
	@cd tests-case-studies; Rscript -e "devtools::load_all(); library(testthat); test_file('test-case-studies_collins-2016.R')"
	@cd tests-case-studies; Rscript -e "devtools::load_all(); library(testthat); test_file('test-case-studies_li-2017.R')"
	@cd tests-case-studies; Rscript -e "devtools::load_all(); library(testthat); test_file('test-case-studies_mcroy-1970.R')"

### ** coverage

## coverage : determine test coverage
.PHONY: coverage
coverage:
	@printf "\n"
	@printf "$(GREEN)=== Determining test coverage ===$(NC)\n"
	@printf "\n"
	@rm -fr docs/coverage/
	@mkdir -p docs/coverage/
	@Rscript -e "library(covr); cov = package_coverage(); report(cov, \"$(TOP_DIR)/docs/coverage/coverage.html\"); cat(paste(\"Coverage_percent: --\", round(percent_coverage(cov), 2), \"--\\n\"))"

### ** pkgdown-offline

## pkgdown-offline : make "pkgdown" (without internet connection)
.PHONY: pkgdown-offline
pkgdown-offline: document update-pkgdown-yml
	@printf "\n"
	@printf "$(GREEN)=== Building the package website with pkgdown (offline) ===$(NC)\n"
	@printf "\n"
	@Rscript -e 'options(pkgdown.internet = FALSE); pkgdown::build_site()'

### ** precompile-doc

## precompile-doc : precompile all Rmd vignettes from Rmd.original files
.PHONY: precompile-doc
precompile-doc:
	@printf "\n"
	@printf "$(GREEN)=== Precompiling vignettes from *.Rmd.original files ===$(NC)\n"
	@printf "\n"
	cd vignettes; Rscript -e 'knitr::knit("tutorial-010-quick-start.Rmd.original", output = "tutorial-010-quick-start.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("tutorial-020-replication.Rmd.original", output = "tutorial-020-replication.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("tutorial-030-steady-state-comps.Rmd.original", output = "tutorial-030-steady-state-comps.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("tutorial-040-pulse-drip-events.Rmd.original", output = "tutorial-040-pulse-drip-events.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("tutorial-050-fixed-effects.Rmd.original", output = "tutorial-050-fixed-effects.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("tutorial-060-units-priors.Rmd.original", output = "tutorial-060-units-priors.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("tutorial-070-prior-predictive-checks.Rmd.original", output = "tutorial-070-prior-predictive-checks.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("tutorial-080-mcmc-output-format.Rmd.original", output = "tutorial-080-mcmc-output-format.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("tutorial-090-post-run-analyses.Rmd.original", output = "tutorial-090-post-run-analyses.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("tutorial-100-posterior-predictive-checks.Rmd.original", output = "tutorial-100-posterior-predictive-checks.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("tutorial-110-derived-parameters.Rmd.original", output = "tutorial-110-derived-parameters.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("tutorial-120-howto-simulations.Rmd.original", output = "tutorial-120-howto-simulations.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("tutorial-130-parameter-identifiability.Rmd.original", output = "tutorial-130-parameter-identifiability.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("case-study-mcroy-1970.Rmd.original", output = "case-study-mcroy-1970.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("case-study-li-2017.Rmd.original", output = "case-study-li-2017.Rmd")'
	cd vignettes; Rscript -e 'knitr::knit("case-study-collins-2016.Rmd.original", output = "case-study-collins-2016.Rmd")'

### ** uncache

## uncache : delete cache files used for vignette rendering
.PHONY: uncache
uncache:
	@printf "\n"
	@printf "$(GREEN)=== Removing cache files ===$(NC)\n"
	@printf "\n"
	@rm -fr vignettes/z-cache-*

### ** check

## check : run 'R CMD CHECK' (also "document")
.PHONY: check
check: document
	@printf "\n"
	@printf "$(GREEN)=== Running 'devtools::check()' ===$(NC)\n"
	@printf "\n"
	@Rscript .run_check_and_get_badge.R

### * --- Package building

## ----PKG BUILDING--------

### ** data

## data : build data files shipped with the pkg
.PHONY: data
data: data/aquarium_mod.rda data/aquarium_run.rda data/eelgrass.rda data/lalaja.rda data/li2017.rda data/li2017_counts.rda data/li2017_prots.rda data/trini_mod.rda
	@printf "\n"
	@printf "$(GREEN)===  Package data files ready ===$(NC)\n"
	@printf "\n"
	@cd prep-data/; rm -f Rplots.pdf

# Models and runs
data/aquarium_mod.rda data/aquarium_run.rda: prep-data/prep-run-example.R
	@printf "\n"
	@printf "$(BLUE)---  Building package data files 'aquarium_*.rda' ---$(NC)\n"
	@printf "\n"
	@mkdir -p data/
	@cd prep-data/; rm -f z-cache-prep-run-example.rds; Rscript prep-run-example.R

data/trini_mod.rda: prep-data/prep-large-model.R prep-data/collins-dump.tsv
	@printf "\n"
	@printf "$(BLUE)---  Building package data file 'trini_mod.rda' ---$(NC)\n"
	@printf "\n"
	@mkdir -p data/
	@cd prep-data/; Rscript prep-large-model.R

# Case studies datasets
data/eelgrass.rda: prep-data/prep-dataset-mcroy1970.R prep-data/mcroy-1970-fig2-topleft.tsv prep-data/mcroy-1970-fig2-topright.tsv prep-data/mcroy-1970-fig2-bottomleft.tsv prep-data/mcroy-1970-fig2-bottomright.tsv
	@printf "\n"
	@printf "$(BLUE)---  Building package data file 'eelgrass.rda' ---$(NC)\n"
	@printf "\n"
	@mkdir -p data/
	@cd prep-data/; Rscript prep-dataset-mcroy1970.R


data/li2017.rda data/li2017_counts.rda data/li2017_prots.rda: prep-data/prep-dataset-li2017.R prep-data/li-2017_design.tsv prep-data/li-2017_proteins.tsv prep-data/li-2017_rel-abundances.tsv prep-data/li-2017_labelling.tsv
	@printf "\n"
	@printf "$(BLUE)---  Building package data files 'li2017*.rda' ---$(NC)\n"
	@printf "\n"
	@mkdir -p data/
	@cd prep-data/; Rscript prep-dataset-li2017.R

data/lalaja.rda: prep-data/prep-dataset-collins2016.R data/trini_mod.rda
	@printf "\n"
	@printf "$(BLUE)---  Building package data file 'lalaja.rda' ---$(NC)\n"
	@printf "\n"
	@mkdir -p data/
	@cd prep-data/; Rscript prep-dataset-collins2016.R

### ** document

## document : generate pkg doc with roxygen2 (also update 'README.md')
.PHONY: document
document: README
	@printf "\n"
	@printf "$(GREEN)=== Generating package documentation with roxygen2 ===$(NC)\n"
	@printf "\n"
	@Rscript -e 'devtools::document(quiet = TRUE)'

### ** [README]

##README : generate README.md from README.Rmd
.PHONY: README
README:
	@printf "\n"
	@printf "$(GREEN)=== Updating 'README.md' from 'README.Rmd' ===$(NC)\n"
	@printf "\n"
	@Rscript -e 'rmarkdown::render("README.Rmd", output_format = rmarkdown::github_document(html_preview = FALSE), quiet = TRUE)'

### ** install

## install : install pkg (also "data", "document", "clean-compiled")
.PHONY: install
install: data document
	@$(MAKE) clean-compiled
	@$(MAKE) onlyinstall
	@printf "\n"
	@printf "$(GREEN)=== Package installed ===$(NC)\n"
	@printf "\n"

### ** [onlyinstall]

##onlyinstall : only install the package (no "data", no "document")
.PHONY: onlyinstall
onlyinstall:
	@printf "\n"
	@printf "$(GREEN)=== Installing the package ===$(NC)\n"
	@printf "\n"
	@Rscript -e 'rstantools::rstan_config()'
	@Rscript -e 'devtools::install(upgrade = FALSE, quick = TRUE, quiet = TRUE)'

### ** pkgdown

## pkgdown : build pkg site (also "document")
.PHONY: pkgdown
pkgdown: document update-pkgdown-yml
	@printf "\n"
	@printf "$(GREEN)=== Building package website with pkgdown ===$(NC)\n"
	@printf "\n"
	@rm -fr docs/*
	@Rscript -e 'pkgdown::build_site()'
	@cp CRAN-version_badge.svg docs/
	@cp dev-version_badge.svg docs/

### ** [update-pkgdown-yml]

##update-pkgdown-yml : generate _pkgdown.yml from _pkgdown.master.yml
.PHONY: update-pkgdown-yml
update-pkgdown-yml:
	@cp pkgdown/_pkgdown.master.yml pkgdown/_pkgdown.yml
	@cd pkgdown/; Rscript .updatePkgdown.R

### ** uninstall

## uninstall : uninstall pkg
.PHONY: uninstall
uninstall:
	@printf "\n"
	@printf "$(GREEN)=== Uninstalling the package ===$(NC)\n"
	@printf "\n"
	@Rscript -e 'tryCatch(remove.packages("isotracer"), error = function(e) {})'

### * --- Cleaning

## ----CLEANING------------

### ** clean

## clean : clean all automatically built files
.PHONY: clean
clean: clean-data clean-man clean-docs clean-vignettes clean-pkgdown-yml clean-compiled
	@rm -f R-CMD-check_badge.svg
	@rm -f R-CMD-check_output.txt
	@rm -f rhub-report_*
	@printf "\n"
	@printf "$(GREEN)=== Cleaned generated files and folders ===$(NC)\n"
	@printf "\n"

### ** clean-compiled

## clean-compiled : delete compiled code
.PHONY: clean-compiled
clean-compiled:
	@printf "$(BLUE)--- Cleaning files related to compiled code ---$(NC)\n"
	@rm -f src/*
	@rm -f R/RcppExports.R

### ** [clean-data]

##clean-data : delete package data generated from prep-data/ scripts
.PHONY: clean-data
clean-data:
	@printf "$(BLUE)--- Cleaning data/ folder ---$(NC)\n"
	@rm -f data/*

### ** [clean-man]

##clean-man : delete documentation generated by Roxygen
.PHONY: clean-man
clean-man:
	@printf "$(BLUE)--- Cleaning man/ folder ---$(NC)\n"
	@rm -f man/*.Rd

### ** [clean-pkgdown-yml]

##clean-pkgdown-yml : delete pkgdown/_pkgdown.yml (generated from master template)
.PHONY: clean-pkgdown-yml
clean-pkgdown-yml:
	@printf "$(BLUE)--- Cleaning pkgdown/_pkgdown.yml file ---$(NC)\n"
	@rm -f pkgdown/_pkgdown.yml

### ** clean-vignettes

## clean-vignettes : clean the vignettes folder
.PHONY: clean-vignettes
clean-vignettes:
	@printf "$(BLUE)--- Cleaning vignettes/ folder ---$(NC)\n"
	@cd vignettes; rm -fr *.Rmd *.html *.pdf *.rds *.RData figures/

### ** clean-docs

## clean-docs : clean pkgdown doc
.PHONY: clean-docs
clean-docs:
	@printf "$(BLUE)--- Cleaning docs/ folder ---$(NC)\n"
	@rm -fr docs/*

### * --- Prerelease check

## ----Prerelease checks----

### ** prerelease-check

# Takes about 37min on a workstation using 4 cores

## prerelease-check : thorough checking before CRAN submission
.PHONY: prerelease-check
prerelease-check:
	@printf "\n"
	@printf "$(GREEN)=== Prerelease check starting ===$(NC)\n"
	@printf "\n"
	@printf "$(BLUE)--- The prerelease check will run: ---$(NC)\n"
	@printf "\n"
	@printf "$(BLUE)      make uninstall $(NC)\n"
	@printf "$(BLUE)      make clean $(NC)\n"
	@printf "$(BLUE)      make install $(NC)\n"
	@printf "$(BLUE)      make test $(NC)\n"
	@printf "$(BLUE)      make precompile-doc $(NC)\n"
	@printf "$(BLUE)      make pkgdown $(NC)\n"
	@printf "$(BLUE)      make coverage $(NC)\n"
	@printf "$(BLUE)      make check $(NC)\n"
	@printf "\n"
	@printf "$(BLUE)--------------------------------------$(NC)\n"
	@printf "\n"
	@$(MAKE) uninstall
	@$(MAKE) clean
	@$(MAKE) install
	@$(MAKE) test
	@$(MAKE) precompile-doc
	@$(MAKE) pkgdown
	@$(MAKE) coverage
	@$(MAKE) check
	@printf "\n"
	@printf "$(GREEN)=== Prerelease check complete ===$(NC)\n"
	@printf "\n"

### ** platforms check

## platforms-check : check on macOS and Windows
.PHONY: platforms-check
platforms-check:
	@printf "\n"
	@printf "$(GREEN)=== Platforms check starting ===$(NC)\n"
	@printf "\n"
	@Rscript -e 'devtools::check_win_release()'
	@Rscript -e 'devtools::check_win_devel()'
	@Rscript -e 'devtools::check_mac_release()'



### * Notes

# To replace the graphic device used in vignettes, use something like (from vignettes/):
# find . -name "*.Rmd" | xargs sed -i 's/fig.align = \"center\", dev = \"jpeg\"/fig.align = \"center\"/g'


### * Attic (2024-05-06: to delete once confident that those are not needed)

# ### ** document

# ## document : generate package documentation using Roxygen (also compiles src files)
# .PHONY: document
# document: quick data README
# 	@printf "\n"
# 	@printf "$(GREEN)*** Generating package documentation with Roxygen ***$(NC)\n"
# 	@printf "\n"
# 	@Rscript -e 'devtools::document(quiet = TRUE)'

# ### ** quick

# ## quick : "quick" install the package (no "data", no "document", does not recompile Cpp sources)
# .PHONY: quick
# quick:
# 	@printf "\n"
# 	@printf "$(GREEN)*** Installing the package ***$(NC)\n"
# 	@printf "\n"
# 	@Rscript -e 'rstantools::rstan_config()'
# 	@Rscript -e 'devtools::install(upgrade = FALSE, quick = TRUE, quiet = TRUE)'

# ### ** buildVignettes

# ## buildVignettes : build the vignettes and put them into inst/doc
# .PHONY: buildVignettes
# buildVignettes: 
# 	@printf "\n"
# 	@printf "$(GREEN)*** Building vignettes with 'devtools::build_vignettes()' ***$(NC)\n"
# 	@printf "\n"
# 	@Rscript -e 'devtools::build_vignettes()'

# ## clean-prep-data-cache : delete cache files in prep_data/
# .PHONY: clean-prep-data-cache
# clean-prep-data-cache:
# 	@rm -f prep-data/z-cache-prep-run-example.rds

# ### ** clean-vignettes-not-cache

# ## clean-vignettes-not-cache : clean the vignettes but not the vignette cache
# .PHONY: clean-vignettes-not-cache
# clean-vignettes-not-cache: 
# 	@rm -fr inst/docs
# 	@cd vignettes; rm -f *.Rmd *.html figures/*

### ** clean-not-cache

# ## clean-not-cache : same as clean but keep vignettes and prep-data cache files
# .PHONY: clean-not-cache
# clean-not-cache: clean-data clean-man clean-docs clean-vignettes-not-cache clean-compiled
# 	@printf "\n"
# 	@printf "$(GREEN)*** Cleaned generated files and folders ***$(NC)\n"
# 	@printf "\n"
# 	@rm -f R-CMD-check_output.txt
# 	@rm -f R-CMD-check_badge.svg

