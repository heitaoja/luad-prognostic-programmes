RSCRIPT ?= Rscript

.PHONY: all check params session
all:
	$(RSCRIPT) --vanilla R/run_pipeline.R

check:
	$(RSCRIPT) --vanilla R/99_verify.R

params session:
	$(RSCRIPT) --vanilla R/00_setup.R
