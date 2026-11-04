# Source-inspection ledger

> **Note added 4 November 2026.** This record dates from the 1 November 2026 search and is kept as written, apart from removing internal file references. Two points now differ from the final manuscript. (1) The final manuscript proves a quantified operator-norm rate (Theorem 2.4: O(ε^{1/2}) in L²₀ operator norm, O(ε) from H² to L²); only the sharp O(ε) operator-norm rate is open. Statements below that no rate is asserted are superseded. (2) The final manuscript credits the operator-level extrinsic coefficient to prior Gaussian/diffusion-operator expansions (including Liou–Lu and, on its abstract, Guo) and claims only the inverse-metric (Sinkhorn metric) consequences. It also adds exact spheres (Proposition 9.1) and embedding dependence (Proposition 10.1). See the addendum at the end of novelty_report.md, which records the 4 November priority check.

**Audit date:** 1 November 2026. **Scope:** the fixed ambient squared-cost theorem and curvature correction in the accompanying manuscript. This ledger records what was actually inspected, rather than treating a bibliography entry as evidence of full-text review.

**Inspection labels.** *Full text, selective* means that the complete paper was available and the specified definitions, statements or proof passages were read; it does not mean that every proof in that source was independently verified. *Team inspection* identifies passages read by another research agent and reported with exact locations. *Abstract only* means that conclusions about the paper's contents remain limited to its abstract.

**Notation alignment.** The manuscript calls the entropic inverse \(B_\varepsilon\), the shape tensor \(Q=\sum_\alpha S_\alpha^2\), and the curvature differential operator \(A_Q=-\rho^{-1}\operatorname{div}(\rho Q\nabla)\). The audited claim is
\[
B_\varepsilon=A^{-1}-\frac{\varepsilon}{4}A^{-1}A_QA^{-1}+O(\varepsilon^2)
\]
on fixed smooth inputs, with parameter derivatives along smooth positive finite-dimensional families. The positive correction uses \(B_\varepsilon A_QB_\varepsilon\). Both \(B_\varepsilon\) and its positive correction converge to \(A^{-1}\) in \(L^2_0\) operator norm, without a claimed general rate. The correction extends to a bounded positive form on all centred \(L^2\) inputs and differs from \(B_\varepsilon\) by \(O(\varepsilon)\) in operator norm. The displayed second-order expansion remains a smooth-input statement; the unit-circle high-frequency example excludes uniform second-order norm accuracy.

## Direct Sinkhorn-geometry sources

### 1. Original metric and formal limit

Lavenant, H., Luckhardt, J., Mordant, G., Schmitzer, B., & Tamanini, L. (2026). The Riemannian geometry of Sinkhorn divergences. *Annales de l'Institut Henri Poincaré C, Analyse non linéaire, 43*(5), 1079–1144. https://doi.org/10.4171/AIHPC/165

**Primary text:** https://arxiv.org/html/2405.04987v2  
**Inspection:** Full text, selective, by root and transport reviewers; journal metadata independently checked. Theorem 3.4 and Section 4 provide the tangent/Hessian representation; Section 4.4, equation (35), gives the formal continuous-density limit. **Role:** established metric identity and motivating unresolved limit, not a contribution of this manuscript.

### 2. Moving finite supports

Mordant, G., & Xu, L. (2026). *Convergence of the Sinkhorn Riemannian metric for finitely supported measures* [Preprint]. arXiv. https://arxiv.org/abs/2609.14651

**Primary text:** https://arxiv.org/pdf/2609.14651  
**Inspection:** Full text, selective, by transport reviewer; abstract independently checked. Introduction and Section 3.3 concern fixed positive atomic weights and smoothly moving distinct atoms, including metric derivatives. **Role:** an earlier positive convergence result in a different base-measure regime; does not establish the present density theorem.

### 3. Computation of geodesic distance

Lavenant, H., Luckhardt, J., & Schmitzer, B. (2026). *Approximation and computation of the geodesic Sinkhorn distance* [Preprint]. arXiv. https://arxiv.org/abs/2610.02007

**Primary text:** https://arxiv.org/html/2610.02007v1  
**Inspection:** Full text, selective; introduction, outline and Remark 1.2 directly checked. Remark 1.2 explicitly fixes \(\varepsilon>0\) and separates the small-regularization question. **Role:** spatial/time discretization and geodesic-computation prior work; no collision with the proposed curvature expansion.

### 4. Hilbert-manifold formulation

Navas, P. (2026). *The information geometry of the Sinkhorn divergence on a Hilbert manifold* [Preprint]. arXiv. https://arxiv.org/abs/2610.39850

**Primary text:** https://arxiv.org/html/2610.39850v1  
**Inspection:** Full text, selective; fixed-regularization geometry, its gap bound, and final discussion inspected. The discussion explicitly leaves general-density metric/Riesz-map/sphere-coordinate limits unresolved. **Role:** recent evidence for the broader open direction, dated 29 October 2026, and a close fixed-parameter geometric framework.

### 5. Second-order geometry and statistical inference

Saidi, T. L., Gunsilius, F., & Mena, G. (2026). *Second-order Sinkhorn geometry and \(\sqrt n\)-inference* [Preprint]. arXiv. https://arxiv.org/abs/2610.36488

**Inspection:** Primary indexed abstract only. It announces differentiability, Levi-Civita and geometry ODEs, and statistical inference. **Role:** nearby fixed-regularization work. Full-text comparison remains incomplete; it is not used as decisive evidence that the target theorem is absent.

## Bistochastic kernels and diffusion maps

### 6. Spectral and semigroup convergence on the torus

Wormell, C. L., & Reich, S. (2021). Spectral convergence of diffusion maps: Improved error bounds and an alternative normalization. *SIAM Journal on Numerical Analysis, 59*(3), 1687–1734. https://doi.org/10.1137/20M1344093

**Primary text:** https://arxiv.org/html/2006.02037v3  
**Inspection:** Full text, selective; torus setting, spectral/semigroup results, and scale analysis inspected. Theorems 3.3 and 3.5 and Proposition 3.6 were additionally checked by the transport reviewer. **Role:** strong prior results for periodised heat-type kernels; precludes novelty claims for the first normalized generator limit.

### 7. Embedded-manifold graph Laplacian

Cheng, X., & Landa, B. (2024). Bi-stochastically normalized graph Laplacian: Convergence to manifold Laplacian and robustness to outlier noise. *Information and Inference: A Journal of the IMA, 13*(4), iaae026. https://doi.org/10.1093/imaiai/iaae026

**Primary text:** https://arxiv.org/html/2206.11386v3  
**Inspection:** Full text, selective; Assumption 1, Lemmas 3.2–3.4, Theorem 4.2, spectral corollary and Section 7 checked. The scale-comparison lemma assumes scale bounds; the main consistency result concerns constrained approximate empirical normalization. **Role:** closest normalization/generator antecedent on embedded manifolds.

### 8. Generator under a prescribed bistochastic measure

Marshall, N. F., & Coifman, R. R. (2019). Manifold learning with bi-stochastic kernels. *IMA Journal of Applied Mathematics, 84*(3), 455–482. https://doi.org/10.1093/imamat/hxy065

**Primary text:** https://arxiv.org/html/1711.06711v2  
**Inspection:** Full text, selective, by root: Theorem 1 and Sections 4.2–4.5 inspected. Journal metadata independently verified. **Role:** foundational generator theorem; full-text access and the assigned proof comparison are no longer gaps. The report makes no claim of novelty for that generator.

### 9. Fixed-bandwidth normalization and sampling refinement

Kessler, N., Magnet, R., & Feydy, J. (2026). Sinkhorn normalization of diffusion kernels. In *Proceedings of the 43rd International Conference on Machine Learning* (Proceedings of Machine Learning Research, Vol. 306, pp. 56764–56795). PMLR. https://proceedings.mlr.press/v306/kessler26a.html

**Primary text:** https://arxiv.org/abs/2507.06161v2  
**Inspection:** Full arXiv v2 PDF, dated 29 May 2026, obtained by the optimisation reviewer and independently read here in Section 2, Section 4/Theorem 4.5, and Appendices B–C. Section 2 explicitly fixes diffusion scale. Appendix B proves symmetric normalization and spectrum in \([0,1]\); Theorem 4.5/Appendix C prove sampling-refinement convergence at fixed bandwidth. **Role:** adjacent normalization/stability result, not an inverse-metric small-bandwidth bias theorem. The earlier abstract-only status is withdrawn. Verified first author: **Nathan**, hence APA initial **N.**

## Gaussian and geometric operator bias

### 10. Ambient Gaussian integral expansion

Liou, J.-M., & Lu, C.-C. (2025). *Asymptotic expansion of the Gaussian integral operators on Riemannian submanifolds of \(\mathbb R^n\)* [Preprint]. arXiv. https://arxiv.org/abs/2506.13238

**Primary text:** https://arxiv.org/html/2506.13238v1  
**Inspection:** Full text, selective; geometric preliminaries, Theorem 3.2, and explicit first curvature coefficient checked. The first author's displayed name is Jia-Ming (Frank) Liou. **Role:** existing Gaussian expansion theory. The existence of a Laplace expansion and raw extrinsic curvature terms is not claimed as new.

### 11. Target-measure diffusion bias

Sule, S., Evans, L., & Cameron, M. (2025). Sharp error estimates for target measure diffusion maps with applications to the committor problem. *Applied and Computational Harmonic Analysis, 79*, 101803. https://doi.org/10.1016/j.acha.2025.101803

**Primary text:** https://arxiv.org/html/2312.14418v1  
**Inspection:** Full text, selective; bias statements, Lemma 3.2, and discussion of the Dirichlet problem checked; journal metadata verified. **Role:** detailed generator bias and inverse-PDE consistency are established nearby results. Its explicit density reweighting differs from exact symmetric Sinkhorn normalization.

### 12. Unresolved full-text comparison

Guo, H. (2026). *Intrinsic and extrinsic bias of diffusion operators: Asymptotic expansion and kernel-dependent comparison* [Preprint]. SSRN. https://doi.org/10.2139/ssrn.7411565

**Primary record:** https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7411565  
**Inspection:** Primary abstract/metadata, posted 4 September 2026, 23 pages. The optimisation reviewer retrieved the landing page; advertised download links returned to that page. **Role:** material unresolved comparison for the raw operator coefficient. The abstract announces intrinsic/extrinsic generator bias at uniform density; it does not settle whether the full paper contains an equivalent logarithmic coefficient or further normalization consequences. This is the main remaining priority uncertainty.

## Heat-kernel metric precedents

### 13. Sharp heat-kernel tensor convergence

Schubr, S. (2026, October 27). *Sharp convergence of heat-kernel Sinkhorn metrics* (Public edition 2) [Author-published preprint]. Luminary Research Archive. https://research.luminaryoddities.com/sharp-convergence-of-heat-kernel-sinkhorn-metrics/v/2

**Primary PDF:** https://research.luminaryoddities.com/sharp-convergence-of-heat-kernel-sinkhorn-metrics/v/2/files/Sharp_Heat_Kernel_Sinkhorn_Metrics.pdf  
**Inspection:** Full text, selective; exact heat-kernel/cost definition, Theorems 2.2–2.3 and principal proof inspected. Theorem 2.2 includes full-\(L^2\) operator-norm convergence at sharp rate \(O(h)\); increasing-frequency examples separate this from improved smooth-input rates. **Role:** prior rigorous smooth-density metric convergence and parameter geometry for a flat-torus, bandwidth-dependent heat cost. It prevents unqualified “first density limit” or “first norm limit” claims.

### 14. Heat-kernel bias and positive correction

Schubr, S. (2026, October 27). *Leading bias and positive correction of heat-kernel Sinkhorn metrics* (Public edition 2) [Author-published preprint]. Luminary Research Archive. https://research.luminaryoddities.com/leading-bias-positive-correction-heat-kernel-sinkhorn-metrics/v/2

**Primary PDF:** https://research.luminaryoddities.com/leading-bias-positive-correction-heat-kernel-sinkhorn-metrics/v/2/files/Leading_Bias_Positive_Correction_Sinkhorn_Metrics.pdf  
**Inspection:** Full text, selective, by optimisation reviewer: Sections 1.1–1.2 and Theorems 1.1–1.2, not all 22 pages of proofs. The model uses \(H_h=e^{h\Delta/2}\), entropy \(2h\), and cost \(c_h=-2h\log[k_h/k_h(0)]\). The theorem has a quadratic leading smooth-input bias with fourth-order remainder; its positive correction is explicit. **Role:** positive correction and entropic-operator sandwiching are prior ideas. The present first-order shape-tensor term belongs to a different fixed-cost model.

## Application and exact-check references

### 15. Natural-gradient interpretation

Shen, Z., Wang, Z., Ribeiro, A., & Hassani, H. (2020). Sinkhorn natural gradient for generative models. In *Advances in Neural Information Processing Systems* (Vol. 33). https://papers.nips.cc/paper_files/paper/2020/hash/122e27d57ae8ecb37f3f1da67abb33cb-Abstract.html

**Primary PDF:** https://papers.nips.cc/paper_files/paper/2020/file/122e27d57ae8ecb37f3f1da67abb33cb-Paper.pdf  
**Inspection:** Full text, selective; Section 3 fixes regularization and Sections 4–5 identify/use the Sinkhorn information matrix. **Role:** application context and earlier Hessian-based optimization, not the present asymptotic theorem.

### 16. Bessel formula and asymptotics

National Institute of Standards and Technology. (n.d.). *Digital Library of Mathematical Functions*, Sections 10.35 and 10.40. Retrieved November 1, 2026, from https://dlmf.nist.gov/10.35 and https://dlmf.nist.gov/10.40

**Inspection:** Equations 10.35.1–10.35.2 and 10.40.1 directly inspected. **Role:** classical generating-function and large-argument identities verifying circle eigenvalues and sharpness constants. No special-function novelty is claimed.

### 17. Discarded quantum collision

Márquez González, S. A. (2026). *Feasibility ordering of entanglement-source placement for qubit channels* [Preprint]. arXiv. https://arxiv.org/abs/2609.18803

**Primary PDF:** https://arxiv.org/pdf/2609.18803  
**Inspection:** Full text, selective; Proposition 1, main feasibility theorem and Section VII inspected. **Role:** decisive exact collision for the abandoned qubit conjecture, including the same transpose-factorization/normal-form proof mechanism. It is retained solely as provenance for the research search.

## Analytic infrastructure added in the norm-limit rerun

### 18. Classical manifold norm-resolvent convergence

Lebeau, G., & Michel, L. (2010). Semi-classical analysis of a random walk on a manifold. *The Annals of Probability, 38*(1), 277–315. https://doi.org/10.1214/09-AOP483

**Primary PDF:** https://arxiv.org/pdf/0802.0644  
**Inspection:** Full text, selective: introduction, geodesic-ball walk and invariant measure definitions, and Theorems 1–2 inspected. Theorem 2 on PDF page 4 gives \(O(h^2)\) \(L^2\) norm-resolvent convergence of the rescaled ball-walk Laplacian in the stated spectral region. **Role:** substantial prior norm-resolvent theory for a different nonlocal kernel and normalization. The present manuscript does not claim a new general norm-resolvent or compactness principle.

### 19. Gaussian heat-kernel comparison

Li, P., & Yau, S.-T. (1986). On the parabolic kernel of the Schrödinger operator. *Acta Mathematica, 156*, 153–201. https://doi.org/10.1007/BF02399203

**Primary reference:** https://doi.org/10.1007/BF02399203  
**Inspection:** Primary full-text Corollary 3.1 selectively inspected by root. **Role:** the heat-kernel Gaussian upper bound used as analytic infrastructure. The weighted operator follows from the manuscript's ground-state transform and bounded-potential comparison. This is not asserted as new.

## Remaining limits

The source ledger does not certify historical priority or independently validate every cited theorem. Guo is the principal unresolved full-text comparison. Saidi et al. remains abstract-level here, with lower apparent relevance to the small-regularization coefficient. The recorded Marshall proof comparison is complete, and the Li–Yau heat-kernel bound was checked in primary full text by root. None of these limitations permits absolute firstness language; the manuscript's claim should remain the precise theorem proved, accompanied by the bounded-search record.
