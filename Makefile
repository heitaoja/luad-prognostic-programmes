RSCRIPT ?= Rscript

.PHONY: all deps fetch check params session

all:
	$(RSCRIPT) --vanilla scripts/run_pipeline.R

deps:
	$(RSCRIPT) --vanilla scripts/install_packages.R

fetch:
	$(RSCRIPT) --vanilla scripts/fetch_inputs.R

check:
	$(RSCRIPT) --vanilla scripts/99_verify.R

params session:
	$(RSCRIPT) --vanilla scripts/00_setup.R
