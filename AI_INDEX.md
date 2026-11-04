# AI index: extrinsic curvature bias in the small-regularisation limit of Sinkhorn metrics (v0.2.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) and the paper
[extrinsic_curvature_sinkhorn_metrics.pdf](extrinsic_curvature_sinkhorn_metrics.pdf) (source [manuscript.tex](manuscript.tex)).
DOI: https://doi.org/10.5281/zenodo.23955055 (concept DOI 10.5281/zenodo.23940147).

## Normalisation

M is a compact connected smooth submanifold of R^D without boundary (dimension m >= 1), with induced metric and volume measure;
mu = rho dx with rho smooth and strictly positive, unit mass. Cost c(x,y) = |x-y|^2 (no factor 1/2), regularisation eps > 0.
S_eps is the Sinkhorn divergence; the Sinkhorn metric g_eps is the limit of S_eps(mu, mu + t h mu)/t^2 (equation (8)), written
g_eps(h mu, h mu) = <h, B_eps h>_mu with B_eps = (eps/2)(I - P_eps^2)^{-1} P_eps on centred L^2(mu) (equation (7)), P_eps the
reversible Markov operator of the exact self-scaling (equations (5)-(6)). A = -rho^{-1} div(rho grad), A phi = h.
Shape operators S_alpha in an orthonormal normal frame; Q = sum_alpha S_alpha^2 (equation (2)); A_Q = -rho^{-1} div(rho Q grad).
On a circle or sphere of radius R, Q = R^{-2}. Exact claims and fingerprints: [prior_art/target.yaml](prior_art/target.yaml).

## Claims and where they are proved

| ID | Claim | Location in the PDF | Evidence type |
|---|---|---|---|
| C1 | Theorem 2.1 (limit and extrinsic bias): (i) B_eps -> A^{-1} in operator norm on centred L^2; (ii) for smooth centred h, B_eps h = A^{-1}h - (eps/4) A^{-1} A_Q A^{-1} h + O(eps^2) in every fixed spatial norm; (iii) g_eps(h mu, h mu) = int \|grad phi\|^2 dmu - (eps/4) int <Q grad phi, grad phi> dmu + O(eps^2). The order eps is sharp on the circle | Section 2 (statement), Section 6 (proof), Section 9.1 (sharpness) | proof; floating-point checks |
| C2 | Corollary 2.2 (positive correction): a positive bilinear form built from B_eps itself, accurate to O(eps^2) on fixed smooth directions, converging to A^{-1} in operator norm | Section 2, proof in Section 6 | proof; floating-point check (ellipse: corrected error at eps = 0.01 is 29.8 times smaller) |
| C3 | Theorem 2.3 (smooth density families): the expansions hold in every fixed parameter-derivative order along smooth positive finite-dimensional density families; the induced finite-dimensional metric G_eps = G_0 - (eps/4)(...) + O(eps^2) | Section 2, proof in Section 8 (Lemma 8.1) | proof |
| C4 | Theorem 2.4 (quantified operator-norm rates): \|\|B_eps - A^{-1}\|\|_{H^2 -> L^2} <= C eps and \|\|B_eps - A^{-1}\|\|_{L^2 -> L^2} <= C eps^{1/2}; same for the corrected operator. Whether O(eps) holds in L^2 operator norm in general is OPEN (Remark 7.2) | Section 2, proof in Section 7 | proof; the sharp rate is open |
| C5 | Lemma 3.1 (relative scaling stability): residual estimate for positive-kernel scaling without pointwise ratio bounds | Section 3 | proof |
| C6 | Lemma 4.1 (gap): uniform spectral gap of I - P_eps of order eps on centred L^2, with uniform bounds | Section 4.3 | proof |
| C7 | Lemma 5.1 (second-order tensor): T = Q/16 in the second-order term of the normalised transition operator; Proposition 5.2: P_eps = I - (eps/4) A + eps^2 (A^2/32 - A_Q/16) + O(eps^3) on smooth functions in every fixed spatial norm (equivalently log P_eps = -eps A/4 - eps^2 A_Q/16 + O(eps^3)). Guo's Theorem 1 with the Gauss equation already gives the tensor Q/16 (implied, not stated, by him; Remark 5.3); credited | Section 5 | proof |
| C8 | Lemma 6.1 (uniform gradient bound) and Proposition 6.2 (operator-norm convergence of the metric operator) | Section 6, Section 6.1 | proof |
| C9 | Lemma 7.1 (uniform Sobolev remainder), the input to Theorem 2.4 | Section 7 | proof |
| C10 | Lemma 8.1 (parameter derivatives and residual upgrading) | Section 8 | proof |
| C11 | Uniform circle: exact eigenvalues b_k = (eps/2) p_k/(1 - p_k^2), p_k = I_k(2R^2/eps)/I_0(2R^2/eps); b_k = (R^2/k^2)(1 - eps/(4R^2) + O(eps^2)) for unit-mass normalisation, i.e. eigenvalue 1/k^2 becomes (1/k^2)(1 - eps/4) on the unit circle | Section 9.1 | proof (Bessel asymptotics, DLMF 10.32, 10.40); floating-point checks |
| C12 | Proposition 9.1 (spheres S^m(R)): the metric computed exactly through the Funk-Hecke formula; explicit eps^2 coefficient ((m-1)^2 + 2L - L^2 - 7)/(96 R^4) of (Lambda b - 1), L = l(l+m-1) | Section 9.2 | proof; symbolic check (code/sympy_sphere_coefficients.py); 60-digit floating-point checks for m = 1..4 |
| C13 | Proposition 9.2 (high-frequency obstruction): on the circle, no L^2 operator-norm rate better than O(eps) for the uncorrected or corrected operator; the eps^2 expansion cannot be uniform in frequency | Section 9.5 | proof; floating-point evidence |
| C14 | Proposition 10.1 (embedding dependence): isometric embeddings with the same density give metrics differing at first order by (eps/4) int <(Q_1 - Q_2) grad phi, grad phi> dmu; Q never vanishes identically (closed manifolds in R^D). Examples: round circle vs ellipse of equal length (coefficients -1/8 and -0.04262 for phi = cos s; averaged -1/8 vs -0.2059); flat tori in R^4 | Section 10 | proof; floating-point checks (code/check_spheres_tori.py) |
| C15 | Examples: product torus in R^4 (anisotropic Q), torus of revolution at nonconstant density, nonuniform ellipse | Sections 9.3, 9.4, 9.6 | floating-point checks only |

Practical reading (Section 10): the relative bias is about eps/(4R^2) for the smallest curvature radius R, so it matters only when
eps is not small compared with 4R^2; intrinsic costs (squared geodesic distance, heat-kernel costs) do not produce this term at
first order; the correction needs the embedding curvature.

## Verification artefacts

Scripts are floating-point checks of selected formulas and constants (the sphere script is a symbolic check); they are not proofs and
not interval certificates. Each exits non-zero on any failed check, also under `python -O` (no assert statements).

- [code/sympy_sphere_coefficients.py](code/sympy_sphere_coefficients.py): symbolic sphere coefficients. Output [results/sympy_sphere_coefficients.json](results/sympy_sphere_coefficients.json).
- [code/check_spheres_tori.py](code/check_spheres_tori.py): spheres S^1-S^4 at 60 digits, Funk-Hecke check, flat torus in R^4, equal-length ellipse vs circle, circle operator-norm evidence. Output [results/check_spheres_tori.json](results/check_spheres_tori.json).
- [code/check_torus_revolution.py](code/check_torus_revolution.py): torus of revolution, nonconstant density. Output [results/check_torus_revolution.json](results/check_torus_revolution.json).
- [code/check_guo_comparison.py](code/check_guo_comparison.py): Guo's Gaussian moments versus the sphere, product-torus, ellipse and torus-of-revolution coefficients, and the first-order identity for every dimension. Output [results/check_guo_comparison.json](results/check_guo_comparison.json).
- [code/check_ellipse_bias.py](code/check_ellipse_bias.py): nonuniform ellipse at 768 and 1536 nodes. Outputs [results/ellipse_bias_results.json](results/ellipse_bias_results.json), [results/ellipse_bias_refinement_768.json](results/ellipse_bias_refinement_768.json), [results/ellipse_bias_refinement_1536.json](results/ellipse_bias_refinement_1536.json).
- [code/independent_circle_check.py](code/independent_circle_check.py), [code/independent_ellipse_check.py](code/independent_ellipse_check.py): independent implementations. Outputs [results/independent_circle_check.json](results/independent_circle_check.json), [results/independent_ellipse_check.json](results/independent_ellipse_check.json).
- [code/make_figures.py](code/make_figures.py): [figures/sinkhorn_curvature_validation.pdf](figures/sinkhorn_curvature_validation.pdf) and [results/figure_validation_summary.json](results/figure_validation_summary.json).
- Outputs reproduce to roundoff, not bit for bit across Python/BLAS builds ([environment.json](environment.json)). Proof notes (not authoritative): [notes/](notes/sinkhorn_manifold_proof.md).

## Exclusions

Population densities only (no empirical-sample or sampling-rate results); closed manifolds without boundary; smooth strictly positive
density; fixed ambient cost |x-y|^2. No convergence of global geodesic distances on the space of measures. The correction needs known
curvature. The sharp O(eps) L^2 operator-norm rate is open on general manifolds.

## Review and provenance

Internal AI-assisted review only: [review/REVIEW_SUMMARY.md](review/REVIEW_SUMMARY.md). Prior art and priority boundary:
[prior_art/novelty_report.md](prior_art/novelty_report.md), [prior_art/source_inspection.md](prior_art/source_inspection.md).
Guo (SSRN 2026, doi:10.2139/ssrn.7411565) was read completely. His Theorem 1 with the Gauss equation gives the tensor Q/16, and his Example 3 with Corollary 1 gives the first sphere and circle bias; the coefficient of the inverse formula follows formally. The contribution here is the rigorous inverse-metric theorem (exact scaling and singular resolvent at arbitrary positive density), not the coefficient. Correction history: [VERSIONS.md](VERSIONS.md). No human peer review, independent reproduction or formal verification. Bibliography: [refs.bib](refs.bib).
