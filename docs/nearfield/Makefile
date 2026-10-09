# Builds the near-field note with Tectonic. Sources in src/, the PDF and the logs in out/.
#
#   make           build out/main.pdf
#   make clean     remove out/
#   make generated regenerate src/generated/*.tex with SymPy (needs a Python with sympy and mpmath, see tools/)

TECTONIC ?= tectonic
PYTHON ?= python3
SRC_DIR := src
OUT_DIR := out
MAIN := main

SOURCES := $(shell find $(SRC_DIR) -type f)

.PHONY: all clean generated

all: $(OUT_DIR)/$(MAIN).pdf

$(OUT_DIR)/$(MAIN).pdf: $(SOURCES) | $(OUT_DIR)
	$(TECTONIC) -X compile $(SRC_DIR)/$(MAIN).tex --outdir $(OUT_DIR) --keep-logs

$(OUT_DIR):
	mkdir -p $@

clean:
	rm -rf $(OUT_DIR)

generated:
	cd tools && $(PYTHON) derive_G3.py && $(PYTHON) make_tex.py && cp G_table.tex explicit_table.tex pieces_tabular.tex ../$(SRC_DIR)/generated/
