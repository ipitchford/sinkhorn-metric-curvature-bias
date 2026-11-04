# Internal review summary

**Paper:** *Extrinsic curvature bias in the small-regularisation limit of Sinkhorn metrics*. Author: Anonymous. Version 0.1.0-candidate.
**Status:** this is an unrefereed, AI-assisted research candidate. The review described here was internal to Evidence Press and was carried out by AI systems. It is not peer review.

## What was done

1. **Five internal reviewer roles.** Five reviewers read the first draft independently, each run as a separate Claude Code subagent: an editor, a methodology (proof-rigour) reviewer, a domain (literature and priority) reviewer, a perspective (significance and readability) reviewer, and a devil's advocate. Each one also ran its own numerical checks, separately from the paper's scripts.
2. **Separate-model adversarial review.** A reviewer from a different model family (OpenAI Codex) audited every result and the prior-art position.
3. **Clean replay.** The lead reran all of the bundle's scripts from a fresh copy.
4. **One revision round**, which included extensions to the results (listed below).
5. **Proof audit.** The new and rewritten proofs were audited line by line. All scripts were rerun from a fresh copy, both normally and with Python's `-O` flag, and the sphere coefficients were checked independently.
6. **Re-review.** Every item on the fix list was checked against the revised manuscript.
7. **Separate-model audit of the new proofs** (OpenAI Codex).

## What was checked

- **Constants and normalisation.** Two independent checks confirmed the Sinkhorn Hessian normalisation: a finite-difference Hessian of a log-domain Sinkhorn solver on a space curve, and the closed-form Bessel eigenvalues on the circle. Several checks confirmed the factor ε/4 and its sign for the cost |x−y|² (with no ½).
- **The key tensor.** Several reviewers re-derived by hand the second-order tensor Q/16 in Lemma 5.1. The derivation covers the Gaussian moments, the chord and volume expansions, the normal-coordinate correction Δ²f = −(2/3)Ric and the Gauss equation. One reviewer also checked it by an exact sum over all 15 Gaussian pairings.
- **Higher-dimensional examples.** On spheres S¹ to S⁴, high-precision von Mises–Fisher eigenvalues reproduce the predicted relative bias −ε/(4R²). This rules out the alternative coefficients that the circle examples alone could not exclude, since on a curve those alternatives coincide with Q. On the flat torus S¹(r₁)×S¹(r₂) in R⁴, the predicted anisotropic bias Q = diag(r₁⁻², r₂⁻²) matches.
- **Every numbered result.** The separate-model review audited each one and judged every one correct under its stated hypotheses. No counterexample was found.

## Substantive objections and how they were handled

- **Priority: Guo (SSRN 2026) and Liou–Lu.** This was the main objection. Guo, *Intrinsic and Extrinsic Bias of Diffusion Operators* (SSRN, doi:10.2139/ssrn.7411565), announces intrinsic and extrinsic bias expansions for radial-kernel diffusion generators at uniform density. Its full text could not be read: SSRN returned HTTP 403, and no other copy was found. It may therefore contain the raw operator-level extrinsic coefficient. Liou–Lu (arXiv 2506.13238) was read in full. It gives a general Gaussian expansion, but its only explicit curvature coefficient is a zeroth-order term in |H|² and scalar curvature.
  - *Disposition:* the paper was repositioned. It now credits the operator-level intrinsic/extrinsic expansions to prior work: Smolyanov–von Weizsäcker–Wittich, Liou–Lu, Guo (on its abstract), Sule–Evans–Cameron, and the classical diffusion-maps papers of Coifman–Lafon, Singer and Hein–Audibert–von Luxburg. It states the structural fact log P_ε = −εA/4 − ε²A_Q/16 + O(ε³). It claims as new only the inverse-metric (Sinkhorn metric) consequences: the formula −(ε/4)A⁻¹A_Q A⁻¹ at arbitrary smooth positive density, with rigorous control of the exact scaling and the singular resolvent, and the results built on it. The inverse-metric statement is described as "not located in the sources we could read". The paper says explicitly that Guo's full text was unread. Historical priority is not asserted.
- **Heat-kernel precedent.** The devil's advocate argued that recent author-posted heat-kernel Sinkhorn papers (Schubr, October 2026) already contain the same overall architecture, with an operator-norm rate.
  - *Disposition:* the paper now says this plainly and contrasts the two settings. It also proves a quantified rate (Theorem 2.4).
- **No operator-norm rate.** The first draft proved operator-norm convergence with no rate.
  - *Disposition:* Theorem 2.4 was added (see below).
- **An unproved sphere statement.** The draft stated the sphere result in a single sentence without proof.
  - *Disposition:* it was replaced by Proposition 9.1, which has a proof and a symbolic certificate.
- **A false example.** One proposed example claimed that a flat torus in R⁴ has Q = 0. Three reviewers showed this is false.
  - *Disposition:* the example was corrected, and the embedding-dependence result was developed into Proposition 10.1.
- **Compressed proofs.** Lemma 5.1 and the parameter-derivative induction in Theorem 2.3 were compressed, and Lemma 5.1 did not state its curvature sign convention.
  - *Disposition:* both proofs were written out and the sign convention was added.
- **Readability and practical use.** Reviewers asked for a plain-English opening, an explanation of what Q means, and practical guidance.
  - *Disposition:* the manuscript gained a plain-English overview ("The paper in brief", with the chord-versus-arc picture and the circle example), plain readings of Q, and practical notes. The practical notes cover the rule of thumb ε ≪ 4R², the point that intrinsic costs avoid this bias, and the fact that the correction needs the curvature to be known.
- **Missing or inexact citations.**
  - *Disposition:* references were added. Every bibliography entry was copied from a Crossref, arXiv or DataCite record.
- **Gaps found by the proof audit and the separate-model audit.** These were: a dropped support restriction and a uniform-injectivity step in Lemma 7.1; lower-order residual terms in the induction for Theorem 2.3; the constant in Proposition 10.1 holding only at a point, not on a neighbourhood; a wording slip about which curve "minimises" the bias; and a stale sentence saying no rate was proved.
  - *Disposition:* all of these were fixed before the final build.

## Extensions added in revision

- **Theorem 2.4 (quantified operator-norm rates).** The error is O(ε) on H² inputs, and O(ε^{1/2}) in L²₀ operator norm for both the Sinkhorn metric operator and its corrected version. The proof uses a new uniform Sobolev remainder estimate (Lemma 7.1).
- **Proposition 9.1 (spheres).** On round spheres S^m(R) with uniform density, the exact eigenvalues are Bessel ratios, and the paper gives an explicit ε² coefficient. This was checked symbolically and at high precision.
- **Proposition 10.1 (embedding dependence).** Two isometric embeddings of the same manifold give Sinkhorn metrics that differ at order ε, and Q is positive definite near a point of maximal radius. The examples are a circle compared with an ellipse of equal length, and flat tori.
- **New check scripts** for spheres, flat tori, the ellipse, and a non-umbilic torus of revolution at non-constant density. The scripts are deterministic and contain no `assert` statements; each one signals failure through its exit code.

## Remaining limits

- The results concern population densities only. There are no sampling or finite-sample rates.
- The manifolds are closed (compact, without boundary) and embedded in Euclidean space. The density must be smooth and strictly positive.
- The positive correction needs the second fundamental form, so the curvature must be known.
- There is no global convergence of geodesic distances. The geometry statements concern smooth finite-dimensional density families.
- The sharp O(ε) operator-norm rate on a general manifold is open. The circle supports it numerically, and the paper's high-frequency obstruction (Proposition 9.2) shows that ε is the best possible order.

## What was not done

- There was no human peer review.
- No unaffiliated party has independently reproduced the work.
- There was no formal (machine-checked) verification of the proofs. The scripts are floating-point and symbolic consistency checks, not proofs.
- The prior-art search was bounded: English-language sources, public indexes, and no subscription databases. No enquiry was sent to any author, including Guo.
