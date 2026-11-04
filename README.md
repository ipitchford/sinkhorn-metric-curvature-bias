# Extrinsic curvature bias in the small-regularisation limit of Sinkhorn metrics

**Unrefereed research candidate, version 0.1.0-candidate (revised 4 November 2026). Author: Anonymous.** Not externally peer reviewed. AI systems assisted the search, derivations, checking, code and writing. No numerical calculation is a premise of a theorem.

DOI: https://doi.org/10.5281/zenodo.23940148 (all versions: https://doi.org/10.5281/zenodo.23940147). Plain-English page: https://evidencepress.org/releases/sinkhorn-metric-curvature-bias/

## In brief

Optimal transport measures how far apart two probability distributions are by the cheapest way of moving one onto the other. Its entropic (Sinkhorn) regularisation, with strength `ε`, is faster to compute, and it also defines a local "metric": how far apart two nearby distributions are. As `ε → 0` this metric should recover the Wasserstein geometry.

For a smooth positive density on a closed manifold embedded in Euclidean space, with fixed cost `|x-y|²`, the paper proves that this metric converges to the Wasserstein metric as `ε → 0` and finds the first error term:

`gε(hμ,hμ) = ∫|∇φ|² dμ − (ε/4)∫〈Q∇φ,∇φ〉 dμ + O(ε²)`,  with `Aφ = h`, `A = −ρ⁻¹div(ρ∇)`, `Q = Σα Sα²`.

The chord between two points is shorter than the arc, and `Q` (squared bending of the manifold) measures exactly this. On a circle or sphere of radius R, Q = R⁻²: on the circle the eigenvalue `1/k²` becomes `(1/k²)(1 − ε/(4R²))` to first order. On the flat product torus S¹(r₁)×S¹(r₂) ⊂ R⁴, Q = diag(r₁⁻², r₂⁻²). Two embeddings of the same intrinsic manifold (a round circle and an ellipse of equal length) give different biases.

Results (numbering as in the PDF):
- operator-norm convergence of the metric operator (Theorem 2.1, Proposition 6.2), with quantified rates O(ε^{1/2}) in L² operator norm and O(ε) from H² to L² (Theorem 2.4); whether the sharp O(ε) operator-norm rate holds in general is open;
- a positive corrected metric accurate to O(ε²) on fixed smooth directions (Corollary 2.2);
- the same expansions along smooth finite-dimensional density families (Theorem 2.3);
- exact round spheres S^m with an explicit ε² coefficient (Proposition 9.1);
- circle high-frequency obstruction: no operator-norm rate better than O(ε) (Proposition 9.2);
- embedding dependence: isometric embeddings give different biases, and Q never vanishes identically (Proposition 10.1).

Credit: the formal limit is due to Lavenant, Luckhardt, Mordant, Schmitzer and Tamanini (2025). Operator-level curvature expansions of Gaussian kernels (Smolyanov–von Weizsäcker–Wittich; Liou–Lu; Guo) are credited. The inverse-metric formula and its consequences were not located in the sources we could read; Guo's SSRN preprint could be read only at abstract level (HTTP 403) and may contain the raw operator-level coefficient. See Section 1 of the paper and `prior_art/`.

## Contents

| File | Purpose |
| --- | --- |
| `extrinsic_curvature_sinkhorn_metrics.pdf` | The paper (23 pages) |
| `manuscript.tex` | LaTeX source (authoritative) |
| `refs.bib`, `refs/` | Bibliography as BibTeX, built by `refs/make_refs.py` from the saved Crossref, DataCite and arXiv records in `refs/records/` |
| `AI_INDEX.md` | Machine-oriented map of claims, evidence and limits |
| `RESEARCH_GATES.json` | Contribution map, searches and limitation triage |
| `review/REVIEW_SUMMARY.md` | Summary of the internal AI-assisted review, objections and their resolution |
| `prior_art/` | Formula-first prior-art search, source-inspection ledger and claim fingerprints (searches of 1 and 4 November 2026) |
| `notes/` | Working proof notes; the manuscript is authoritative |
| `code/`, `results/` | Check scripts and their JSON outputs |
| `figures/` | Validation figure |
| `environment.json`, `requirements.txt` | Replay environment |
| `MANIFEST.sha256` | SHA-256 of every file |
| `LICENSE.md`, `CITATION.cff`, `.zenodo.json`, `.github/workflows/replay.yml` | Licences, citation metadata, archive metadata, CI |

## Reproduce

Python 3 with the pinned packages in `requirements.txt`. All scripts are deterministic and write JSON to `results/`; each exits with a non-zero status on failure (no `assert`, so they also fail correctly under `python -O`). Run them on a copy if you want to keep the shipped outputs.

```bash
export OPENBLAS_NUM_THREADS=1
U="uv run --with numpy==2.3.5 --with scipy==1.17.0 --with matplotlib==3.10.8 --with mpmath --with sympy python3 -I -B"
$U code/sympy_sphere_coefficients.py   # symbolic check of the sphere coefficients
$U code/check_spheres_tori.py          # spheres S^1..S^4 (60 digits), Funk-Hecke check, flat torus in R^4, equal-length ellipse vs circle, circle operator-norm evidence
$U code/check_torus_revolution.py      # torus of revolution, nonconstant density, modes n=0,1
$U code/check_ellipse_bias.py          # nonuniform-ellipse experiment (768/1536 nodes)
$U code/independent_circle_check.py; $U code/independent_ellipse_check.py; $U code/make_figures.py
```

Each script runs in a few seconds. Outputs reproduce to floating-point roundoff, not bit for bit across Python and BLAS builds (see `environment.json`).

Rebuild the PDF: `pdflatex manuscript.tex` twice.

## Numerical checkpoints

- Ellipse (axes 1.7, 0.8): g₀ = 0.48506748107430087, predicted coefficient −0.32545935958807615; corrected metric error at ε = 0.01 is 29.8 times smaller than raw.
- Spheres: ε² coefficient of (Λ b − 1) matches ((m−1)² + 2L − L² − 7)/(96R⁴) to ≥ 8 digits for m = 1..4, ℓ = 1..3.
- Flat torus r = (1, 0.6): slopes match −Λ_Q/(4Λ²) to < 10⁻⁶.
- Torus of revolution: Richardson slopes −0.1601009 / −0.0391023 vs predicted −0.1601019 / −0.0391025.
- Circle: sup_k |b_k − 1/k²|/ε ≈ 0.2500 for ε = 10⁻²…10⁻⁴.

These are floating-point checks, not certified interval bounds.

## Scope and limits

Closed, connected, compact, smoothly embedded manifold; smooth strictly positive density; population densities (no empirical sampling rates). The correction needs the embedding curvature to be known. No results on boundaries, vanishing densities, noisy samples or global geodesic distances.

## Status and provenance

Internal AI-assisted review only (five reviewer roles, a separate-model adversarial review, a clean replay, one revision round, a proof audit, a re-review and a separate-model audit of the new proofs; see `review/REVIEW_SUMMARY.md`). No human peer review, independent reproduction or formal verification. The publication replay reran all seven scripts with and without `python -O` from fresh copies.

## How to cite

Anonymous (2026). *Extrinsic curvature bias in the small-regularisation limit of Sinkhorn metrics*. Version 0.1.0-candidate. Evidence Press. https://doi.org/10.5281/zenodo.23940148

## Licence

Text, figures, data and records CC0-1.0; code MIT. See `LICENSE.md`.
