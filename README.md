# Aging-like BAC sectors when coupling travels through barrier-limited desmoplastic edges

**Thesis #32.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Depends on:** Thesis #25 (BAC sectors versus Gompertz gain) and Thesis #21 (metastasis graph barrier conductances).

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

Do aging-like Bounded Adaptive Coherence sectors that appear under lumped coupling still appear when the same subsystem couples only through barrier-limited desmoplastic edges, or does the barrier schedule carve a different sector map?

The barrier schedule carves a different map when the stalled fractions differ by sector. The shared object is the five-index weight matrix of Thesis #25, read by that deposit's sector classifier. Under the barrier reading each edge is a series pair in the sense of Thesis #21, and the grounded Laplacian sees only the series flux. At a uniform stalled fraction, global decay stays aging-like, because the flux is a common factor and the sector ratio remains 1. Assigning the transport-limited fraction 0.705882 to one sector and the shedding-limited fraction 0.264706 to the other classes that same path neither (sector ratios 0.432781 and 2.310638). On the 169-node rate plane the aging-like sets are then disjoint: 11 lumped nodes and 19 barrier nodes, intersection empty. Co-decay of conductance with shedding restores the lumped map on all 169 nodes. The organism-scale cut stays cancer-like.

No number is taken from either deposit's results file. This deposit does not rescore a Gompertz residual and does not claim rejuvenation or cure.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/barrier_sectors.py` | Shared toy: lumped and barrier-coupled sector maps (seed 20260921) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Eigenvalue paths, rate planes, disagreement, barrier schedule |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/barrier_sectors.py
python3 build_pdf.py
```

NumPy and Matplotlib are required for the toy. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Aging-like BAC sectors when coupling travels through barrier-limited desmoplastic edges [Internet]. Thesis #32 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-32-bac-sectors-under-barrier-transport

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
