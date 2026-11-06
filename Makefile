# bulk / TCGA -- reproducible workflow
# ----------------------------------------------------------------
# R is invoked with an explicit interpreter so the workflow also runs on
# Windows where `Rscript` is not always on PATH.  Override with
#   make RSCRIPT=/path/to/Rscript all
#
# Steps are thin scripts; each writes its tables to results/tables/ and is
# safe to rerun on its own.

RSCRIPT ?= Rscript

STEPS := 00_setup 01_load_data 02_score 03_tcga_primary 03b_mutual_adjustment \
         04_external_meta 05_incremental 06_negative_control 07_figures 08_tables

.PHONY: all clean dirs params session check $(STEPS)

all: $(STEPS)

dirs:
	@mkdir -p data/interim data/processed results/objects results/tables results/figures results/logs

params:
	$(RSCRIPT) -e "source('R/00_setup.R'); write_params()"

session:
	$(RSCRIPT) -e "writeLines(capture.output(sessionInfo()), 'docs/sessionInfo.txt')"

$(STEPS): | dirs
	$(RSCRIPT) R/$@.R 2>&1 | tee results/logs/$@.log

# sanity checks before every commit: no empty tables, params recorded
check:
	@$(RSCRIPT) -e "f <- list.files('results/tables', '\\\\.csv$$', full.names=TRUE); \
	  z <- f[file.size(f) < 60]; if (length(z)) { cat('SUSPICIOUS (tiny) tables:\n'); print(z); quit(status=1) } else cat('tables ok:', length(f), '\n')"

clean:
	rm -rf results/objects results/figures results/tables results/logs data/interim data/processed
