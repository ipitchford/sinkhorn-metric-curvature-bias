# Novelty-gate report: fixed-cost Sinkhorn metric and extrinsic curvature bias

> **Note added 4 November 2026.** This record dates from the 1 November 2026 search and is kept as written, apart from removing internal file references. Two points now differ from the final manuscript. (1) The final manuscript proves a quantified operator-norm rate (Theorem 2.4: O(ε^{1/2}) in L²₀ operator norm, O(ε) from H² to L²); only the sharp O(ε) operator-norm rate is open. Statements below that no rate is asserted are superseded. (2) The final manuscript credits the operator-level extrinsic coefficient to prior Gaussian/diffusion-operator expansions (including Liou–Lu and, on its abstract, Guo) and claims only the inverse-metric (Sinkhorn metric) consequences. It also adds exact spheres (Proposition 9.1) and embedding dependence (Proposition 10.1). See the addendum at the end of novelty_report.md.

**Search date:** 1 November 2026.  
**Researcher/agent:** an AI research agent performing an independent literature and algebra audit for the research team.  
**Corpus boundary:** English-language primary papers retrievable through two live web-search indexes, arXiv HTML/PDF, journal and proceedings pages, authors' public manuscript repositories, SSRN metadata, NIST DLMF, and OEIS. Searches covered material indexed through the search date. No subscription MathSciNet/zbMATH/Scopus search or author enquiry was undertaken.  
**Intended contribution units:** a rigorous theorem for exact population Sinkhorn normalization with fixed ambient squared cost; its explicit first curvature correction; a positive correction with second-order smooth-direction accuracy; and an independently checkable proof/verification package. Historical priority is not asserted.

## 1. Frozen claim

Let \(M^d\subset\mathbb R^D\) be a compact connected smooth embedded manifold without boundary, with its induced metric. Let \(\mu=\rho\,d\mathrm{vol}\), where \(\rho>0\) is smooth and has integral one. The ground cost is the **fixed** function \(c(x,y)=|x-y|^2\), and the entropy multiplier is \(\varepsilon\):

\[
\mathrm{OT}_\varepsilon(\mu,\nu)
=\inf_{\pi\in\Pi(\mu,\nu)}\left\{\int c\,d\pi+\varepsilon\,\mathrm{KL}(\pi\mid\mu\otimes\nu)\right\}.
\]

Let \(u_\varepsilon>0\) satisfy

\[
u_\varepsilon(x)\int_M e^{-|x-y|^2/\varepsilon}u_\varepsilon(y)\,d\mu(y)=1,
\]

and let \(P_\varepsilon\) be the resulting reversible Markov operator. On \(L^2_0(\mu)\), put

\[
B_\varepsilon=\frac\varepsilon2P_\varepsilon(I-P_\varepsilon^2)^{-1},\qquad
A=-\rho^{-1}\operatorname{div}(\rho\nabla),\qquad G=A^{-1}.
\]

The already established Sinkhorn Hessian identity identifies the tangent quadratic form with \(\langle h,B_\varepsilon h\rangle_\mu\). For a local orthonormal normal frame, write \(S_\alpha\) for the shape operators and set

\[
Q=\sum_\alpha S_\alpha^2,\qquad
A_Q=-\rho^{-1}\operatorname{div}(\rho Q\nabla).
\]

The principal candidate statement, matching Theorem 2.1 and the positive-correction corollary of the accompanying manuscript, is, for each smooth mean-zero \(h\),

\[
B_\varepsilon h=Gh-\frac\varepsilon4GA_QGh+O(\varepsilon^2),
\]

and hence, with \(\phi=Gh\),

\[
g_\varepsilon(h\mu,h\mu)
=\int_M|\nabla\phi|^2\,d\mu
-\frac\varepsilon4\int_M\sum_\alpha|S_\alpha\nabla\phi|^2\,d\mu
+O(\varepsilon^2).
\]

The manuscript establishes the remainder in each prescribed smooth norm for fixed smooth inputs, uniformly with every prescribed finite number of parameter derivatives in compact smooth positive density families. It additionally proves the qualitative operator-norm limit
\[
\|B_\varepsilon-G\|_{L^2_0(\mu)\to L^2_0(\mu)}\longrightarrow0.
\]
This includes strong convergence for every centred \(L^2\) direction. No quantitative general rate is asserted for this operator-norm limit. The leading smooth-input coefficient is nonpositive in the quadratic form. It is not claimed to be strictly negative for all nonzero directions on every embedding.

The positive correction, defined on the full centred \(L^2\) space through its bounded bilinear form, is

\[
\widehat B_\varepsilon
=B_\varepsilon+\frac\varepsilon4B_\varepsilon A_Q B_\varepsilon.
\]

For each fixed \(\varepsilon>0\), its quadratic form is strictly positive on nonzero centred \(L^2\) inputs. It satisfies
\[
\|\widehat B_\varepsilon-B_\varepsilon\|_{L^2_0\to L^2_0}\le C\varepsilon,
\qquad
\|\widehat B_\varepsilon-G\|_{L^2_0\to L^2_0}\longrightarrow0.
\]
For each fixed smooth input its error is \(O(\varepsilon^2)\). This is an application of the explicit bias formula, not a claim that the general idea of positive bias correction is new. The uniform gradient estimate used to extend the correction to all \(L^2\), and the compact-domination argument for norm convergence, use standard covariance and spectral methods; neither general method is claimed as new.

### Exclusions

- No bounded Euclidean domain with boundary, boundary layer, vanishing density, noncompact support, or singular measure is covered.
- No general quantitative rate for either full \(L^2\to L^2\) operator-norm limit is claimed. A uniform second-order rate is excluded by the circle example below.
- No unrestricted infinite-dimensional curvature convergence or global geodesic-distance convergence is claimed. Geometry statements concern smooth finite-dimensional families.
- The limiting tangent metric is the intrinsic weighted Wasserstein metric on \(M\). This does not identify ambient-cost transport between arbitrary distant measures with intrinsic geodesic transport.
- The underlying Sinkhorn Hessian formula, diffusion-map generator limit, Gaussian Laplace expansion, Bessel diagonalisation, and positivity of Gaussian kernels are prior mathematics.

## 2. Normalisation and translation table

| Source/object | Canonical notation | Checked translation |
|---|---|---|
| Gibbs kernel for \(c=\lVert x-y\rVert^2\), entropy \(\varepsilon\) | \(e^{-\lVert x-y\rVert^2/\varepsilon}\) | Local Gaussian covariance is \(\varepsilon I/2\); first generator term is \(\varepsilon\Delta/4\). |
| Normalized ambient Gaussian integral | \(T_\varepsilon=(\pi\varepsilon)^{-d/2}\int e^{-c/\varepsilon}(\cdot)\,d\mathrm{vol}\) | \(v_\varepsilon=(\pi\varepsilon)^{d/4}u_\varepsilon\), so \(v_\varepsilon T_\varepsilon(v_\varepsilon\rho)=1\). |
| Heat-kernel manuscripts using \(H_h=e^{h\Delta/2}\) | \(h=\varepsilon/2\) | Their tensor \(hP_h(I-P_h^2)^{-1}\) has the same scalar normalization, but their cost depends on \(h\). |
| Liou–Lu kernel \(e^{-\lVert x-y\rVert^2/(4t)}(4\pi t)^{-d/2}\) | \(t=\varepsilon/4\) | Coefficients must be rescaled by \(4^{-j}\) at order \(j\). |
| Weighted backward diffusion generator | \(-A=\Delta+\nabla\log\rho\cdot\nabla\) | \(P_\varepsilon=I-\varepsilon A/4+O(\varepsilon^2)\). |
| Shape convention | \(\langle S_\alpha X,Y\rangle=\langle\mathrm{II}(X,Y),n_\alpha\rangle\) | \(\sum S_\alpha^2\) is independent of normal-frame rotations and the sign convention for each normal. |
| Mean-curvature notation in coefficient calculations | \(H=\sum_i\mathrm{II}(e_i,e_i)\) | This is the unaveraged trace. Papers defining the average mean curvature require a factor \(d\). |
| Markov Green operator | \(P(I-P^2)^{-1}=\sum_{m\ge0}P^{2m+1}\) on \(L^2_0\) | It is an odd-step potential kernel; the inverse metric is a rescaled discrete Poisson operator. |

## 3. Formula fingerprints

### 3.1 Uniform embedded circle

For a circle of radius (R), the exact normalized transition kernel is the von Mises convolution kernel, with concentration \(\kappa=2R^2/\varepsilon\). Its Fourier eigenvalues are

\[
p_k=I_k(\kappa)/I_0(\kappa),\qquad
r_{\varepsilon,k}=\frac{\varepsilon p_k}{2(1-p_k^2)},\quad k\ne0.
\]

The Bessel generating function is the classical Jacobi–Anger identity

\[
e^{\kappa(z+z^{-1})/2}=\sum_{k\in\mathbb Z}I_k(\kappa)z^k.
\]

This is not a rational or algebraic generating function; canonical rational-denominator searches therefore do not apply. The named special-function fingerprint was searched instead. NIST DLMF §§10.35 and 10.40 confirm the identity and fixed-order large-argument asymptotics.

For \(R=1\), the expansions are

\[
p_k=1-\frac{k^2}{4}\varepsilon+
\frac{k^2(k^2-2)}{32}\varepsilon^2+O(\varepsilon^3),\qquad
r_{\varepsilon,k}=\frac1{k^2}-\frac\varepsilon{4k^2}+O(\varepsilon^2).
\]

For \(k=1,\ldots,6\), the integer fingerprints are \(k^2=(1,4,9,16,25,36)\) and \(k^2(k^2-2)=(-1,8,63,224,575,1224)\). The metric leading values are (1,1/4,1/9,1/16,1/25,1/36); the linear coefficients are (-1/4,-1/16,-1/36,-1/64,-1/100,-1/144). These are classical Laplacian/Bessel fingerprints, not new sequences. The squares match OEIS A000290. The mode \(k=1\) already establishes that the proposed \(O(\varepsilon)\) uncorrected rate cannot generally improve to \(O(\varepsilon^2)\).

The norm-limit strengthening prompted a second fingerprint at increasing frequency. On the unit circle, choose \(\varepsilon_n=4/n^2\) and \(k=n\). Rescaling the exact von Mises integral gives \(p_n\to e^{-1}\), hence
\[
\frac{r_{\varepsilon_n,n}}{\varepsilon_n}\longrightarrow\frac1{4\sinh1}.
\]
Here \(G\) has eigenvalue \(1/n^2=\varepsilon_n/4\) and \(A_Q=A\). The correction adds \(\varepsilon_n n^2r_{\varepsilon_n,n}^2/4=r_{\varepsilon_n,n}^2=O(\varepsilon_n^2)\). Therefore
\[
\frac{\widehat r_{\varepsilon_n,n}-1/n^2}{\varepsilon_n}
\longrightarrow\frac1{4\sinh1}-\frac14<0.
\]
The corrected norm error is thus at least \(c\varepsilon_n\) on this sequence. This does not contradict second-order accuracy on fixed smooth inputs. High-frequency obstructions to improved full-\(L^2\) rates already occur in the heat-kernel literature; the circle computation is an exact verification for this fixed-cost model, not a new general obstruction principle.

### 3.2 Small exact kernel instances

The theorem concerns a continuum operator, not a finite matrix family. Nonetheless the following exact finite uniform-polygon kernels fix indexing and normalization. They are algebra checks, not continuum approximations with a proved error bound.

- Two antipodal unit-circle points: with \(a=e^{-4/\varepsilon}\), \(P=(1+a)^{-1}\begin{pmatrix}1&a\\a&1\end{pmatrix}\), eigenvalues (1,(1-a)/(1+a)).
- Three equally spaced points: with \(b=e^{-3/\varepsilon}\), \(P=[(1-b)I+b\mathbf1\mathbf1^T]/(1+2b)\), eigenvalues (1,(1-b)/(1+2b),(1-b)/(1+2b)).
- Four equally spaced points: with \(t=e^{-2/\varepsilon}\), (P) is circulant with first row \((1,t,t^2,t)/(1+t)^2\), eigenvalues \(1,q,q,q^2\), where \(q=(1-t)/(1+t)\).

All have unit row/column sums, positive entries and positive spectrum for finite positive bandwidth. The finite objects are not claimed new and were not used to establish historical priority of the continuum curvature formula.

### 3.3 Geometric coefficient fingerprint

For \(T_\varepsilon=I+\varepsilon Q_1+\varepsilon^2Q_2+\cdots\), the diagnostic coefficient is the principal second-order part of \(Q_2-Q_1^2/2\):

\[
\frac1{16}\sum_\alpha S_\alpha^2.
\]

Independent Wick contraction reproduces it. Before converting ordinary normal-coordinate fourth derivatives to covariant \(\Delta^2\), the anisotropic Hessian term is

\[
-\mathrm{Ric}/24+\sum_\alpha[(\operatorname{tr}S_\alpha)S_\alpha/48+S_\alpha^2/24].
\]

The conversion adds \(\mathrm{Ric}/48\). Gauss' equation \(\mathrm{Ric}=\sum_\alpha[(\operatorname{tr}S_\alpha)S_\alpha-S_\alpha^2]\) leaves \(\sum S_\alpha^2/16\). Exact self-adjoint normalization then removes the scalar potential and produces the weighted divergence form. This is the formula-level fingerprint around which the final rerun was conducted.

## 4. Cross-domain alias map

| Domain | Aliases tested or explicitly compared |
|---|---|
| Analysis | Nonlocal Laplacian, inverse graph Laplacian, resolvent, Green operator, weighted Poisson equation, semiclassical Gaussian integral expansion. |
| Probability | Reversible Markov kernel, bistochastic transition, odd-step potential kernel, diffusion approximation, ground-state/Doob transform. |
| Harmonic analysis and representation theory | Fourier convolution on the circle, von Mises kernel, modified Bessel ratios, Jacobi–Anger series, circle Laplacian spectrum. |
| Combinatorics | Weighted walk sums \(\sum P^{2m+1}\), square sequence and transformed coefficient sequences. No nontrivial new enumerative object is claimed. |
| Mathematical physics and geometry | Short-time heat kernel, Schrödinger kernel, extrinsic/intrinsic bias, second fundamental form, Weingarten/shape operator, logarithm of a transfer operator. |
| Optimisation, statistics and ML | Sinkhorn Hessian, Sinkhorn information matrix, natural gradient, diffusion maps, target-measure diffusion maps, kernel normalization, curvature bias correction. |

## 5. Reproducible search log

All queries below were executed on 1 November 2026. The two live indexes are denoted S1 and S2. S1 was used for high-recall collision checks; S2 for routine metadata and named-object verification. Negative searches mean that no exact collision was returned and verified, not that no such paper exists.

| Ring | Corpus | Exact query or retrieval | Consequence |
|---|---|---|---|
| Original and citing literature | Primary arXiv/journals | Lavenant et al. `2405.04987`; Mordant–Xu `2609.14651`; Navas `2610.39850`; Saidi et al. `2610.36488` | Density limit is identified as formal/open in the first, second and Navas sources; latest geometry paper concerns fixed regularization. |
| Adjacent diffusion maps | S1/S2 and full primary HTML | `"Spectral convergence of diffusion maps"`; `"Bi-stochastically normalized graph Laplacian"`; `"Manifold learning with bi-stochastic kernels"` | Generator/spectral convergence is established prior art and cannot constitute the novelty claim alone. |
| Exact new object | S1 | `Sinkhorn metric extrinsic curvature shape operator Gaussian kernel second order expansion` | Located Gaussian expansion and latest geometric papers, no exact inverse-metric curvature formula. |
| Exact terminology | S1 | `"Sinkhorn" "second fundamental form"`; `"Sinkhorn" "extrinsic" "bias"`; `"Sinkhorn" "shape operators"`; `"Sinkhorn" "second fundamental" metric bias`; `"Sinkhorn" "extrinsic curvature"` | No verified exact collision; heat-kernel bias manuscript was prominent and was inspected. |
| Broader geometry aliases | S1 | `"entropic" "metric tensor" "extrinsic"`; `"Sinkhorn metric" "embedding"`; `"Gaussian" "logarithm" "shape operator" "manifold" integral kernel` | No verified exact collision. |
| Alternative normalization literature | S1 | `"bistochastic" "second fundamental form"`; `bistochastic Gaussian kernel "second-order" manifold curvature`; `"bistochastic" "curvature" "second order"`; `"bistochastic" "logarithm" "curvature" Gaussian` | General Gaussian bias work and normalization literature located; none inspected states the target Hessian formula. |
| New adjacent bias literature | S1/S2 | `diffusion maps second order bias extrinsic curvature shape operators kernel Gaussian`; `"Intrinsic and Extrinsic Bias of Diffusion Operators"`; `"Hongfu Guo" "diffusion" "bias"`; `"Hongfu Guo" "second fundamental" "kernel"` | Located Guo SSRN7411565, a material abstract-only gap described below. |
| Gaussian operator coefficients | S1 and primary HTML | `Liou Lu 2506.13238 Gaussian integral operator embedded manifold second order`; primary `arxiv.org/html/2506.13238v1`; primary `arxiv.org/html/2312.14418v1` | Existing general expansion and target-diffusion bias must be credited; neither inspected source states the Sinkhorn inverse bias. |
| Markov/Green aliases | S2/S1 | `"Sinkhorn" "resolvent" "Wasserstein"`; `"Sinkhorn metric" "Green" operator`; `"bistochastic" "inverse" "Laplace" "metric"`; `"Sinkhorn" "Dirichlet form" "curvature"` | Navas is a close fixed-parameter source; no small-noise shape-operator result located. |
| ML/natural-gradient aliases | S1, NeurIPS PDF, PMLR | `"Sinkhorn Natural Gradient" "Wasserstein" convergence`; `"Sinkhorn Normalization of Diffusion Kernels" curvature asymptotic` | Shen et al. use fixed regularization; Kessler et al. full text confirms fixed-bandwidth normalization and sampling refinement (§4, Appendices B–C). |
| Harmonic fingerprint | S2, DLMF | `"von Mises" "I_n" "I_0" convolution eigenvalues`; `"Sinkhorn" "I_k" "I_0"`; direct DLMF10.35 and10.40 | Bessel diagonalisation is classical. No associated Sinkhorn curvature theorem located. |
| Sequence fingerprint | S2/OEIS | `"1, 4, 9, 16, 25, 36" site:oeis.org`; direct OEIS A000290 | Confirms classical square sequence; gives no independent novelty evidence. |
| Positive-correction rerun | Primary author archive | `leading-bias-positive-correction-heat-kernel-sinkhorn-metrics/v/2` | Earlier positive corrections exist for a different heat-kernel model. Claim only the explicit fixed-cost correction derived from the shape tensor. |

The operator-norm strengthening triggered the following additional searches, also on 1 November 2026:

| Ring | Corpus | Exact queries | Consequence |
|---|---|---|---|
| Norm/resolvent structural rerun | S1 | "Sinkhorn metric" "norm" "resolvent" manifold; "Sinkhorn" "operator norm" convergence metric manifold; "bistochastic" "norm resolvent" convergence; "diffusion maps" "norm resolvent" compact manifold | Reidentified Schubr's prior heat-model norm theorem and classical nonlocal norm-resolvent theory. No exact fixed-cost metric collision was returned and verified. |
| Classical norm-resolvent comparison | S2 and primary arXiv PDF | "Lebeau" "Michel" "2010" "random" compact manifold spectrum; "nonlocal Laplacian" "norm resolvent" "manifold"; "Sinkhorn" "inverse weighted Laplacian" "norm"; "bistochastic" "inverse" "operator norm"; "Semi-classical analysis of a random walk on a manifold" arxiv | Full selective inspection of Lebeau–Michel, Theorems 1–2, confirms quantitative norm-resolvent convergence is established for geodesic-ball walks; different kernel and normalization. |
| High-frequency and compactness aliases | S2 | "Sinkhorn" "sinh" "operator norm"; "Sinkhorn" "high frequency" "correction"; "Sinkhorn" "compact domination" | No further exact fixed-cost theorem was returned and verified. Existing heat-model sharpness remains prior art. |
| Corrected full-\(L^2\) extension rerun | S1 | "Sinkhorn" "conditional variance" gradient operator; "Sinkhorn" "Fisher" "operator norm" correction; "bistochastic" "gradient" "resolvent" norm | Reidentified the heat-model positive-correction source. The row-score covariance estimate is treated as standard analysis, not a new method. |

The exact queries are sufficient to rerun the bounded search; search rankings and indexed versions can change. Source passages below, rather than a query's failure, carry the positive evidence for scope distinctions.

## 6. Verified prior results and closest comparisons

### 6.1 Originating metric and later convergence papers

**Hugo Lavenant, Jonas Luckhardt, Gilles Mordant, Bernhard Schmitzer, Luca Tamanini.** *The Riemannian geometry of Sinkhorn divergences*. Annales de l'Institut Henri Poincaré C 43(5) (2026), 1079–1144. DOI: <https://doi.org/10.4171/AIHPC/165>. Preprint: <https://arxiv.org/abs/2405.04987>. Primary journal page: <https://ems.press/journals/aihpc/articles/14299213>. Journal metadata inspected here; the transport agent inspected the full PDF, especially §4.4, equation (35). The Hessian/Markov inverse formula is their result. Their smooth-density small-regularization discussion is formal and motivates the current target. No priority claim is made for the metric itself.

**Gilles Mordant, Liane Xu.** *Convergence of the Sinkhorn Riemannian metric for finitely supported measures*, arXiv:2609.14651, 13 September 2026. <https://arxiv.org/abs/2609.14651>. Abstract inspected here; full PDF inspected by the transport agent. Their regime has fixed positive atomic weights and smoothly moving distinct atoms. It includes metric-derivative and curvature conclusions in that finite-support setting. It does not prove the continuous-density theorem considered here. This is important prior evidence, not a claim that every later extension remains open.

**Paul Navas.** *The information geometry of the Sinkhorn divergence on a Hilbert manifold*, arXiv:2610.39850v1, 29 October 2026. Full primary HTML inspected: <https://arxiv.org/html/2610.39850v1>. The final discussion, near rendered line692, explicitly says that general-density metric/Riesz-map/sphere-coordinate convergence remains to be proved. This is the most recent inspected explicit statement of the broader gap. The paper treats fixed-regularization geometry on positive Sobolev density strata. Its fixed-parameter gap bound does not replace the \(c\varepsilon\) gap needed here.


**Hugo Lavenant, Jonas Luckhardt, Bernhard Schmitzer.** *Approximation and computation of the geodesic Sinkhorn distance*, arXiv:2610.02007v1, 1 October2026. Full primary HTML inspected: <https://arxiv.org/html/2610.02007v1>. Remark1.2 explicitly fixes the regularization parameter and separates its spatial/time discretization analysis from the small-regularization question. It is therefore relevant geometric-computation prior art, not a collision with the present limit.

### 6.2 The closest diffusion-map results

**Caroline L. Wormell, Sebastian Reich.** *Spectral Convergence of Diffusion Maps: Improved Error Bounds and an Alternative Normalization*. SIAM Journal on Numerical Analysis 59(3) (2021), 1687–1734. DOI: <https://doi.org/10.1137/20M1344093>. Full primary HTML inspected: <https://arxiv.org/html/2006.02037v3>. The setting is a flat hypertorus with a periodised Gaussian kernel. Spectral and semigroup convergence, including improved bias under Sinkhorn normalization, strongly precede qualitative inverse convergence in that model. Their operator results are central prior art; a bare inverse-limit theorem in their exact setting would be weak novelty. A periodised heat kernel on a flat torus has no ambient shape tensor of the target kind.

**Xiuyuan Cheng, Boris Landa.** *Bi-stochastically normalized graph Laplacian: convergence to manifold Laplacian and robustness to outlier noise*. Information and Inference 13(4) (2024), iaae026. DOI: <https://doi.org/10.1093/imaiai/iaae026>. Full primary HTML inspected: <https://arxiv.org/html/2206.11386v3>; journal metadata verified at <https://academic.oup.com/imaiai/article/13/4/iaae026/7762892>. Their assumptions include a compact connected embedded manifold without boundary and positive smooth sampling density. Lemma3.3 compares approximate Sinkhorn scales using PSD positivity but assumes pointwise scale bounds. Theorem4.2 gives empirical generator consistency for a constrained approximate scaling procedure. Their spectral corollary assumes crude low-spectrum convergence; §7 identifies stronger scale consistency as a route to further spectral analysis. The present candidate's unconditional residual estimate and exact population inverse asymptotics should be distinguished from these results, not presented as the first generator limit.

**Nicholas F. Marshall, Ronald R. Coifman.** *Manifold learning with bi-stochastic kernels*. IMA Journal of Applied Mathematics 84(3) (2019), 455–482. DOI: <https://doi.org/10.1093/imamat/hxy065>. <https://arxiv.org/abs/1711.06711>. Full primary HTML was opened by the root reviewer: <https://arxiv.org/html/1711.06711v2>. Theorem1 and §§4.2–4.5 were inspected; its infinitesimal generator under prescribed bistochastic measures is established prior art. That proof comparison is complete at the recorded selective-inspection level. No inverse-metric shape-bias theorem was identified in the inspected material.

**Gilles Lebeau, Laurent Michel.** *Semi-classical analysis of a random walk on a manifold*. Annals of Probability 38(1) (2010), 277–315. DOI: <https://doi.org/10.1214/09-AOP483>. Full primary PDF selectively inspected: <https://arxiv.org/pdf/0802.0644>, introduction and Theorems1–2, especially Theorem2 on PDF page4. For the averaging walk on intrinsic geodesic balls of radius \(h\), their rescaled discrete Laplacian has \(O(h^2)\) \(L^2\) norm-resolvent convergence in the stated spectral region. The invariant measure involves ball volume and is not a prescribed smooth density under Gaussian Sinkhorn normalization. This is substantial classical adjacent work: norm-resolvent convergence, and finite-spectral-projection proofs of compact limits more generally, are not new principles.

**Peter Li, Shing-Tung Yau.** *On the parabolic kernel of the Schrödinger operator*. Acta Mathematica156 (1986), 153–201. DOI: <https://doi.org/10.1007/BF02399203>. The root reviewer inspected the primary Corollary3.1 used for the Gaussian heat-kernel upper bound. The weighted comparison in the manuscript follows after the usual ground-state transform and bounded-potential comparison. This source is analytic infrastructure, not a novelty competitor.

### 6.3 Existing Gaussian/extrinsic operator expansions

**Jia-Ming (Frank) Liou, Chi-Chien Lu.** *Asymptotic Expansion of the Gaussian Integral Operators on Riemannian submanifolds of \(\mathbb R^n\)*, arXiv:2506.13238v1, 16 June 2025. Full primary HTML inspected: <https://arxiv.org/html/2506.13238v1>. Theorem3.2 gives a general expansion; the first coefficient is explicitly expressed through mean and scalar curvature. The paper credits earlier ambient-versus-geodesic distance expansions to Smolyanov, von Weizsäcker and Wittich. It does not treat Sinkhorn normalization or inverse metrics. Finite-order smooth Laplace asymptotics used in the current proof are standard; no novelty claim should be attached to their existence.

**Shashank Sule, Luke Evans, Maria Cameron.** *Sharp error estimates for target measure diffusion maps with applications to the committor problem*. Applied and Computational Harmonic Analysis 79 (2025), 101803. DOI: <https://doi.org/10.1016/j.acha.2025.101803>. Full primary preprint HTML inspected: <https://arxiv.org/html/2312.14418v1>. The article derives explicit bias prefactors for target-measure diffusion maps and a second-order kernel expansion (Lemma3.2), then controls a Dirichlet BVP. Its normalization is an explicit density reweighting, not the exact symmetric Sinkhorn scale. It establishes that detailed generator-bias and inverse-PDE consistency work already exists in adjacent literature.

**Hongfu Guo.** *Intrinsic and Extrinsic Bias of Diffusion Operators: Asymptotic Expansion and Kernel-Dependent Comparison*. SSRN7411565, posted 4 September 2026, 23 pages. DOI: <https://doi.org/10.2139/ssrn.7411565>. Primary metadata/abstract: <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7411565>. The indexed primary abstract was inspected; full text could not be retrieved (403). It announces explicit intrinsic/extrinsic generator bias for radial kernels on closed submanifolds with uniform density, including circle/sphere comparisons. This is a material unresolved comparison for the operator-level \(\sum S_\alpha^2\) coefficient. The abstract does not announce exact Sinkhorn normalization, an inverse metric, arbitrary positive density, or the positive correction above. Absence from the abstract does not prove absence from the paper.

### 6.4 Very close recent heat-kernel results

**Scott Schubr.** *Sharp convergence of heat-kernel Sinkhorn metrics*. Luminary Research Archive, public edition2, 27 October 2026; manuscript dated21 October2026. Primary page: <https://research.luminaryoddities.com/sharp-convergence-of-heat-kernel-sinkhorn-metrics/v/2>. Full seven-page PDF inspected: <https://research.luminaryoddities.com/sharp-convergence-of-heat-kernel-sinkhorn-metrics/v/2/files/Sharp_Heat_Kernel_Sinkhorn_Metrics.pdf>. Theorem2.2 supplies sharp norm/Sobolev tensor estimates on a flat torus for \(H_h=e^{h\Delta/2}\); Theorem2.3 treats parameter derivatives and finite-family geometry. The cost is \(c_h=-2h\log[k_h/k_h(0)]\), with entropy (2h). It explicitly separates that model from a fixed Euclidean-cost problem. This is publicly posted prior work irrespective of its unrefereed status. It precludes claiming the first rigorous smooth-density Sinkhorn metric limit without the fixed-cost qualification.

**Scott Schubr.** *Leading bias and positive correction of heat-kernel Sinkhorn metrics*. Same archive, public edition2, 27 October2026. Primary page and full PDF inspected by the optimisation reviewer, especially §§1.1–1.2 and Theorems1.1–1.2: <https://research.luminaryoddities.com/leading-bias-positive-correction-heat-kernel-sinkhorn-metrics/v/2/files/Leading_Bias_Positive_Correction_Sinkhorn_Metrics.pdf>. It announces a second-order heat-model bias involving the Wasserstein Hessian of Fisher information and a positive correction with fourth-order smooth-family accuracy. This makes “a positive correction of Sinkhorn geometry” an existing idea. The proposed contribution here is the different first-order, always nonpositive, extrinsic-shape term for a fixed ambient Gaussian cost, and its associated \(O(\varepsilon^2)\) correction. It is therefore treated as inspected prior work, not as an abstract-only gap; the comparison concerns the different kernel and bias order.

The norm-limit rerun explicitly compared Theorem2.2 of the first Schubr paper: its full-\(L^2\) operator-norm rate is \(O(h)\), with sharp high-frequency examples. The present qualitative operator-norm limit is therefore scoped to the fixed-cost embedded-manifold model. Neither the existence of a norm limit nor the loss of smooth-input accuracy at increasing frequency is claimed as a general first result.

### 6.5 Other immediate aliases

**Zebang Shen, Zhenfu Wang, Alejandro Ribeiro, Hamed Hassani.** *Sinkhorn Natural Gradient for Generative Models*. NeurIPS33 (2020). Full primary conference PDF inspected: <https://papers.nips.cc/paper_files/paper/2020/file/122e27d57ae8ecb37f3f1da67abb33cb-Paper.pdf>. The Sinkhorn information matrix is the divergence Hessian and is used for natural gradients. Section3 explicitly keeps the regularization parameter fixed. This supports application relevance but not the target small-regularization asymptotics.

**Nathan Kessler, Robin Magnet, Jean Feydy.** *Sinkhorn Normalization of Diffusion Kernels*. ICML2026, PMLR306, 56764–56795. Primary proceedings metadata: <https://proceedings.mlr.press/v306/kessler26a.html>. The full arXiv2507.06161v2 PDF (29 May2026) was obtained by the optimisation reviewer and its relevant passages independently inspected here: <https://arxiv.org/abs/2507.06161v2>. Section2 explicitly fixes bandwidth; Theorem4.5 studies weak sampling refinement at fixed scale. AppendixB proves symmetric normalization and spectrum in [0,1]; AppendixC proves the continuous scaling and sampling limit. These are normalization/stability results, not a small-regularization inverse-metric bias theorem. The earlier access gap is closed.

**Tristan Luca Saidi, Florian Gunsilius, Gonzalo Mena.** *Second-order Sinkhorn Geometry and \(\sqrt n\)-Inference*, arXiv:2610.36488, 27 October2026. Primary indexed abstract inspected: <https://arxiv.org/abs/2610.36488>. It concerns metric differentiability, connections, geometry ODEs and statistical inference. Full text was unavailable through the retrieval tool. No small-regularization assertion was found in the abstract.

**NIST Digital Library of Mathematical Functions.** §§10.35,10.40, <https://dlmf.nist.gov/10.35>, <https://dlmf.nist.gov/10.40>. Both pages inspected. They fix the classical Bessel generating-function and asymptotic identities used in the circle check. **OEIS A000290**, <https://oeis.org/A000290>, inspected for the square-spectrum fingerprint.

### 6.6 Proof mechanism behind the norm-limit strengthening

The manuscript's heat-kernel comparison gives, as quadratic forms on the centred space,
\[
I-P_\varepsilon\ge c(I-e^{-c_0\varepsilon A}).
\]
Positivity and inverse order yield
\[
0\le B_\varepsilon
=\frac{\varepsilon}{4}\bigl[(I-P_\varepsilon)^{-1}-(I+P_\varepsilon)^{-1}\bigr]
\le C(G+\varepsilon I).
\]
Let \(F_N\) project onto the first \(N\) positive eigenmodes of \(A\), with eigenvalues \(\lambda_j\), and let \(K_N=I-F_N\). Then the high-frequency block has norm at most \(C(\lambda_{N+1}^{-1}+\varepsilon)\), while positivity bounds the mixed block by its square root times a uniform constant. Strong convergence handles the finite block. Taking \(\varepsilon\to0\), then \(N\to\infty\), proves norm convergence without a rate. These are familiar compactness and operator-order arguments; the model-specific scale and heat comparison make them applicable here.

For the corrected tensor, centering the spatial derivative of a kernel row gives a covariance identity. Its row-score second moment is at most \(C/\varepsilon\), so conditional Cauchy–Schwarz and reversibility give
\[
\|\nabla P_\varepsilon f\|_2^2
\le \frac C\varepsilon
\langle f,(I-P_\varepsilon^2)f\rangle_\mu.
\]
Apply this to \(f=(I-P_\varepsilon^2)^{-1}h\), and use the \(c\varepsilon\) spectral gap, to obtain \(\|\nabla B_\varepsilon h\|_2\le C\|h\|_2\). Thus \(B_\varepsilon A_QB_\varepsilon\) extends to a bounded positive form uniformly in small \(\varepsilon\). Its contribution to \(\widehat B_\varepsilon\) is \(O(\varepsilon)\) in norm. The covariance identity, Gaussian moment bound and spectral-gap argument are standard tools, not separately advertised contributions.

## 7. Discarded quantum conjecture: decisive collision

The earlier candidate was the qubit source-placement implication

\[
\Psi\circ\Phi\text{ non-entanglement-breaking}
\implies \Phi\otimes\Psi\text{ non-entanglement-annihilating}.
\]

The team independently found a proof using a qubit Sinkhorn/Pauli normal form, an invertible completely-positive transpose factorization, and depolarizing approximation of the boundary. Formula-level and high-recall search then found an exact earlier solution:

**Samuel A. Márquez González.** *Feasibility Ordering of Entanglement-Source Placement for Qubit Channels*, arXiv:2609.18803, 16 September2026. Full primary six-page PDF inspected: <https://arxiv.org/pdf/2609.18803>. Proposition1 supplies the same type of transpose factorization; the main theorem proves the implication. Its SectionVII leaves a stronger quantitative SDP conjecture open. The candidate theorem/proof novelty therefore **FAILED** and was abandoned. An all-intermediate-filter equivalence obtained by the same machinery is too close a routine extension to justify a release on its own.

The original source was **Masajada, Fellous-Asiani and Streltsov**, *Optimizing entanglement distribution via noisy quantum channels*, Physical Review A113 (2026), 052414; DOI <https://doi.org/10.1103/hdzn-fwpj>; full preprint HTML inspected at <https://arxiv.org/html/2506.06089v1>. Its conjecture label was not evidence that the statement remained open in November2026. The discarded proof and collision warning were retained in the internal research record.

## 8. Missingness and priority limits

The central uncertainty is **historical priority**, not a counterexample found to the proposed theorem. The explicit continuous-density gap is documented by very recent sources, but that does not imply that every special case or the desired curvature coefficient is unprecedented.

The main material gap is Guo's full operator-bias manuscript. Saidi et al.'s latest full text was also unavailable to this agent, although its indexed abstract concerns fixed-regularization geometry and statistics. The earlier Kessler, Marshall and Schubr-bias access gaps were subsequently closed by the team, with the inspected portions recorded above and in `source_inspection.md`. Older heat-kernel/quantum-path-integral literature may contain an equivalent \(\sum S_\alpha^2\) logarithmic coefficient. A publication should therefore credit Gaussian/extrinsic expansion work and identify the **normalized inverse/Hessian** formula as the proposed novelty, rather than claiming discovery of the raw coefficient. No author enquiry was sent, as external messaging was not authorized or requested. This omission is recorded rather than treated as an obstacle to continued mathematical work.

No finite search can certify historical firstness. Search results in English can omit dissertations, non-English sources, unindexed preprints, changed versions and unpublished arguments. The very recent quantum collision demonstrates why this report avoids absolute novelty language.

## 9. Outcome by contribution unit

| Unit | Outcome | Assessment | Permitted wording |
|---|---|---|---|
| Sinkhorn Hessian/Markov inverse identity | COLLISION / established prior result | Lavenant et al.; earlier information-matrix work | “Using the known Hessian representation…” |
| First diffusion generator \(P_\varepsilon=I-\varepsilon A/4+\cdots\) | COLLISION / established prior result | Bistochastic diffusion-map literature | “We give the scale control needed to use this generator expansion at the inverse level.” |
| Smooth density limit for heat-kernel flat torus | COLLISION | Wormell–Reich operator results and Schubr's explicit tensor result | Distinguish this model by its bandwidth-dependent heat cost. |
| Exact fixed-cost closed-manifold inverse limit in \(L^2\) operator norm, with full scale/gap control | CLEAR WITHIN RECORDED BOUNDARY; bounded priority uncertainty | No identical inspected theorem found; heat-model norm limits and general norm-resolvent methods are prior | “We prove qualitative operator-norm convergence in this fixed-cost, smooth positive density setting…” |
| Shape-operator first bias for the normalized inverse metric | BOUNDED UNCERTAINTY; strongest candidate unit | No exact match found; raw-operator overlap with Guo remains uninspected | “The resulting expansion identifies the coefficient as …”; “No matching inverse-metric formula was found in the search described here.” |
| Positive correction with smooth-input \(O(\varepsilon^2)\) accuracy and qualitative \(L^2\) norm convergence | CLEAR WITHIN RECORDED BOUNDARY as a specific corollary | Positive correction as a general idea is already prior art; boundedness uses a standard row-score variance estimate | “This curvature coefficient yields the positive corrected tensor …” |
| General compact-domination, heat comparison and covariance machinery | Established analytic methods; no separate novelty claim | Lebeau–Michel and heat-kernel theory are major antecedents; the present proof adapts familiar tools | Credit sources and state the model-specific application without claiming a new norm-resolvent principle. |
| All smooth finite-family parameter derivatives/curvature | BOUNDED UNCERTAINTY | Similar conclusions already exist for heat kernels and moving atoms; the fixed-cost density scope is the contribution | State the precise density-family theorem without claiming the first curvature convergence of Sinkhorn geometry. |
| Exact circle check and reproducibility artifacts | Verification/exposition contribution | Bessel machinery is classical | “An independent circle check confirms constants and sharpness.” |
| Qubit feasibility conjecture | COLLISION | Exact September2026 theorem and same proof architecture | Cite González; make no new-solution claim. |

Forbidden language includes “the first rigorous density limit”, “the general open problem is solved”, “the most important open problem in applied mathematics”, “no previous work studied this bias”.

## 10. Decision and rerun triggers

**Decision: proceed with the sharply scoped fixed-cost theorem and explicit extrinsic inverse-metric expansion, under BOUNDED UNCERTAINTY about priority.** The current evidence supports a serious research candidate, not a certified historical-first result. Its value is materially stronger than a qualitative generator-to-inverse corollary because it identifies the sign, tensor and sharp order of the embedding-dependent bias and gives a positive correction.

Before submission, the most useful additional search is obtaining Guo's full paper and checking whether the logarithmic Gaussian coefficient is already isolated there; then finish any remaining detailed comparison with the already retrieved normalization and heat-kernel correction papers. Rerun the formula gate if the proof produces a simpler exact spectral function, a known nonlocal Dirichlet form, a heat-kernel conjugacy, or a correction identical to a published operator preconditioner. None of those future checks is replaced by the negative searches recorded above.

## Addendum, 4 November 2026: priority check and public-discussion search

This addendum records two checks made after the internal review, on 4 November 2026. The sections above are left as written on 1 November.

**Guo (2026), read at abstract level only.** Hongfu Guo, *Intrinsic and Extrinsic Bias of Diffusion Operators: Asymptotic Expansion and Kernel-Dependent Comparison*, SSRN, doi:10.2139/ssrn.7411565. The Crossref metadata record (indexed 4 September 2026, 22 references) was saved and its abstract read in full. The full text could not be obtained: the SSRN abstract and download routes returned HTTP 403, and no arXiv or repository copy was found by title, author or keyword search. No attempt was made to bypass the access restriction.

What the abstract announces:

- the continuous bias of the generator of a radial-kernel diffusion operator on a closed Riemannian submanifold, at **uniform** sampling density, with a kernel expansion to order ε⁴;
- a split of the leading bias into intrinsic parts (radial moments, curvature) and extrinsic parts (kernel moments, the second fundamental form and its derivatives), and the discrepancy between geodesic and ambient distance;
- a moment-based criterion for comparing kernel profiles, simplifications for constant curvature and spheres, and circle and sphere numerics.

What the abstract does not announce: exact Sinkhorn (bistochastic, self-scaling) normalisation; the Sinkhorn metric or any inverse (resolvent) of the operator; arbitrary positive density; operator-norm convergence of the metric operator; a positive corrected metric; derivative estimates along density families; or a high-frequency obstruction. Absence from an abstract does not show absence from the paper.

**Decision.** The raw operator-level extrinsic coefficient (a Σ_α S_α²-type term in a Gaussian generator expansion) is plausibly anticipated by Guo and by older Gaussian and heat-kernel expansions (Liou–Lu, Smolyanov–von Weizsäcker–Wittich, and the diffusion-maps literature). The final manuscript therefore credits operator-level intrinsic/extrinsic bias expansions to that prior work. It claims only the inverse-metric consequences: the Sinkhorn metric formula −(ε/4)A⁻¹A_Q A⁻¹ at arbitrary smooth positive density, with control of the exact scaling and of the singular resolvent, together with the results built on it. Its wording is that the inverse-metric statement was "not located in the sources we could read", and it states explicitly that Guo's full text could not be read and may contain the raw coefficient.

**Public-discussion and corpus search.** GitHub issues and discussions were searched for the same objects (Sinkhorn metric, entropic metric, extrinsic curvature, shape operator, second fundamental form). There were no issue hits, and the discussion hits were unrelated (for example optimal-transport library usage questions). No public prior statement of the result was found. The Evidence Press release feed (239 items) was scanned by keyword. No release concerns the Sinkhorn metric tensor or its small-regularisation limit, so this release does not duplicate earlier work. The closest release, *poincare-sinkhorn-global-rates* (doi:10.5281/zenodo.23929617), studies global convergence rates of the Sinkhorn algorithm. That is a different question with no shared theorem, and the manuscript cites it in one sentence.

**Still not done.** No author enquiry was sent, and no subscription database (MathSciNet, zbMATH, Scopus) was searched.
