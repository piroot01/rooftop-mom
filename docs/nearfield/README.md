# Near-field singularity extraction (technical note)

A manual to the near-field extraction in `src`: the three-term extraction (1/R, R, R^3) on all touching cells of the
rooftop EFIE on a uniform grid, and the closed forms of all static integrals it needs, in the form of the
bachelor-thesis appendix (in `a` and `b`, one integer prefactor, `\phantom`-aligned continuation lines).

- `make` builds `out/main.pdf` (Tectonic); `make clean` removes `out/`.
- `make generated` regenerates `src/generated/*.tex` (needs a Python with sympy and mpmath, e.g.
  `make generated PYTHON="uv run --no-project --with sympy --with mpmath python"`).
- `src/include/template.sty` is the bachelor-thesis template (unchanged); `src/include/note.sty` adds only macros.
  Layout (`article`, A4, 1.5in margins) and colours (`Paired*`, `elmagGay`) are those of the thesis.
- `src/figures/`: the TikZ figures (cells and rooftops, wedge pairs, the correlation substitution, the integration
  regions, the corner difference); `\noteroofx` in `note.sty` draws a rooftop.
- `tools/`: `derive_G3.py` (the mixed antiderivatives, verified symbolically and by 30-digit quadrature; writes
  `G3.pkl`, `G3_table.tex` and the MATLAB `nearFieldG.m`), `doc_pieces.py` (pieces, corner sums, symbolic proof of the
  agreement with the thesis formulas), `make_tex.py` (Appendices A and B and the piece table, each formula rebuilt
  from its printed parts and checked to 30 digits). The files they write in `tools/` are not tracked.
