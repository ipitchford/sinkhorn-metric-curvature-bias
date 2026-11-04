# Fixed-cost Sinkhorn metric limit on closed embedded manifolds

Research proof notes, 1 November 2026. These are original derivations in this session; literature checking is separate. The main open target is the rigorous continuous-density limit proposed in Lavenant–Luckhardt–Mordant–Schmitzer–Tamanini (2025/2026), §4.4, equation (35). Mordant–Xu (September 2026) settle moving finite supports and explicitly leave densities at a formal level.

## 1. Statement and conventions

Let M be a compact, connected, C-infinity manifold of intrinsic dimension d, without boundary, smoothly embedded in R^D with induced Riemannian metric. Write dx for its volume measure. Let rho>0 be C-infinity with integral rho dx=1 and mu=rho dx. The ground cost is the FIXED ambient squared distance c(x,y)=|x-y|^2. It is not an epsilon-dependent heat-kernel cost.

Define OT_epsilon(mu,nu)=inf_pi int c dpi + epsilon KL(pi | mu tensor nu), and the debiased Sinkhorn divergence S_epsilon. For bounded mean-zero h, the affine path mu_t=(1+t h)mu is positive for sufficiently small |t|, and its tangent quadratic form is g_epsilon_mu(h mu,h mu)=lim_{t->0} S_epsilon(mu,mu_t)/t^2. For arbitrary L2 mean-zero h, use the continuous extension of this quadratic form supplied by the operator formula below; an unbounded h need not define a positive two-sided affine path.

Let u_epsilon be the unique positive self-Sinkhorn scaling satisfying

    u_epsilon(x) int exp(-c(x,y)/epsilon) u_epsilon(y) dmu(y)=1.

The associated reversible Markov operator is

    P_epsilon f(x)=u_epsilon(x) int exp(-c(x,y)/epsilon) u_epsilon(y) f(y) dmu(y).

The Gaussian kernel restricted to M is positive definite, so P_epsilon is positive semidefinite, self-adjoint in L2(mu), Markov, and a contraction. Constants are its simple eigenvectors at eigenvalue one, because its integral kernel is everywhere positive.

For mean-zero h in L2(mu), the established Hessian formula is

    g_epsilon_mu(h mu,h mu)
      = (epsilon/2) <h, (I-P_epsilon^2)^(-1) P_epsilon h>_mu.

Let A=-Delta_M-grad(log rho).grad be the nonnegative weighted Laplacian, with inverse on mean-zero functions. Define

    R_epsilon=(epsilon/2)(I-P_epsilon^2)^(-1)P_epsilon.

**Theorem A (continuous-density limit and curvature bias).** For every mean-zero h in L2(mu),

    R_epsilon h -> A^(-1)h in L2(mu),

and hence g_epsilon_mu(h mu,h mu)-><h,A^(-1)h>_mu. For smooth h, with phi=A^(-1)h and Q=sum_alpha S_alpha^2 defined from the embedding's shape operators,

    g_epsilon_mu(h mu,h mu)
      = int |grad phi|^2 dmu
        -(epsilon/4) int <Q grad phi,grad phi> dmu
        +O(epsilon^2).

The resolvent applied to smooth h has a full asymptotic expansion in every spatial C^k norm. Constants are uniform over compact smooth families of strictly positive densities and smooth directions.

**Theorem B (finite-dimensional geometry).** For a C-infinity family rho_theta>0 on a finite-dimensional parameter domain, the matrix

    G_epsilon,ij(theta)=g_epsilon_{mu_theta}(partial_i mu_theta,partial_j mu_theta)

converges in C^k on compact parameter sets to the intrinsic Wasserstein pullback metric G_0 for every finite k, with O(epsilon) error. If the density map is an immersion, G_0 is positive definite; consequently the inverse metrics, Levi-Civita symbols, and curvature tensors converge in every corresponding C^k norm. In fact these quantities possess asymptotic power series in epsilon.

The limit is the intrinsic Wasserstein tangent metric on M. Ambient squared-distance optimal transport between distant measures on M is not generally the intrinsic Wasserstein distance; no global-distance identification follows from this theorem.

## 2. A positive-kernel scaling-stability lemma

This is the key elementary estimate that avoids assuming uniform derivatives of the unknown Sinkhorn scaling when using Laplace expansion.

Let (X,mu) be a probability space, K a symmetric positive-semidefinite integral operator with strictly positive kernel, and u,q positive functions. Assume uKu=1 pointwise. Let r=qKq, delta=||r-1||_infinity<1, s=u/q, and Q=diag(q) K diag(q). Then Q1=r and Qs=1/s. The operator Q remains positive semidefinite.

Set a=int s dmu and b=int 1/s dmu. By symmetry,

    b=<1,Qs>=<s,Q1>=int sr dmu.

The Cauchy–Schwarz inequality for the semidefinite bilinear form induced by Q gives

    b^2 <= <s,Qs><1,Q1> = int r dmu,

because <s,Qs>=int 1 dmu=1. Since b>= (1-delta)a,

    a <= sqrt(1+delta)/(1-delta).

Let E=int (s-1)^2/s dmu. Positivity gives

    0 <= <s-1,Q(s-1)>
       = int (s-1)(1/s-r) dmu
       = -E + int (s-1)(1-r) dmu.

Therefore

    E <= delta int |s-1| dmu
      <= delta sqrt(E a).

If E=0 there is nothing to prove; otherwise divide by sqrt(E). Thus

    int (s-1)^2/s dmu
      <= delta^2 sqrt(1+delta)/(1-delta).

No pointwise ratio bound between u and q is assumed. The result is equally valid for finite matrices with positive weights.

## 3. Upgrading stability directly to uniform derivative estimates (root improvement)

The L2-to-uniform interpolation argument below can be replaced by a considerably cleaner argument. Let q be a uniformly positive and bounded approximate *normalised* scaling, let Q=diag(q) T_epsilon diag(q) act on L2(mu), and let s=v_epsilon/q. Then Qs=1/s, Q1=r. The stability lemma implies

    ||s-1||_L1 <= sqrt(E int s) <= C delta.

Because ||Q||_{L1(mu)->C^k}<=C_k epsilon^(-d/2-k) (differentiate only its explicit, smooth kernel and q),

    ||1/s-r||_{C^k}=||Q(s-1)||_{C^k}
       <= C_k epsilon^(-d/2-k) delta.

If r=1+O(epsilon^N) in C^k, choosing N>d/2+k ensures 1/s is uniformly close to one and its derivatives are controlled. Composition by the reciprocal then gives

    ||s-1||_{C^k} <= C_k epsilon^(N-d/2-k).

Thus arbitrary-order approximate scalings give arbitrary-order uniform spatial derivative expansions WITHOUT any prior pointwise/derivative bounds on the unknown exact scaling. This is the preferred final proof. The polynomial estimates below remain useful for parameter-derivative bookkeeping but are unnecessary for spatial asymptotics.

### Optional independent polynomial bounds

Normalise the Gaussian kernel by

    T_epsilon f(x)=(pi epsilon)^(-d/2) int_M exp(-c(x,y)/epsilon) f(y) dy,
    v_epsilon=(pi epsilon)^(d/4) u_epsilon.

Then v_epsilon T_epsilon(v_epsilon rho)=1. Compactness and Gaussian localisation give a uniform bound on T_epsilon 1. Since derivatives of c in x are uniformly bounded, differentiation of the positive integral shows

    |grad log v_epsilon| <= C/epsilon.

If v_epsilon achieves maximum V at x_0, then v_epsilon(y)>=V exp(-C) on the geodesic ball B_M(x_0,epsilon). The ambient distance is at most the geodesic distance, so the Gaussian there is at least exp(-epsilon). A uniform lower bound on rho and the small-ball volume give

    T_epsilon(v_epsilon rho)(x_0) >= c V epsilon^(d/2).

The scaling equation gives 1>=c V^2 epsilon^(d/2), hence

    sup v_epsilon <= C epsilon^(-d/4).

At any x, T_epsilon(v_epsilon rho)(x)<=C sup v_epsilon, whence

    inf v_epsilon >= c epsilon^(d/4).

For every integer j, repeated differentiation of v=1/T_epsilon(v rho) gives

    ||v_epsilon||_{C^j} <= C_j epsilon^(-d/4-j).

To justify the last statement, derivatives of the Gaussian kernel divided by the Gaussian are bounded by C_j epsilon^(-j), because all derivatives of c are bounded. Thus the derivative of the positive integral relative to the integral has the same bound. Faà di Bruno's formula for its reciprocal then gives the stated estimate. Only polynomial bounds matter; these exponents are deliberately crude.

## 4. Arbitrary-order approximate scalings

In local normal coordinates, Laplace expansion gives, for each N,

    T_epsilon f = sum_{j=0}^N epsilon^j L_j f + O(epsilon^(N+1)),

with smooth differential operators L_j, uniformly in any fixed C^k norm provided enough derivatives of f are bounded. Here L_0=I and

    L_1 f=(1/4)Delta_M f + E_M f,

where E_M is a smooth curvature/extrinsic-embedding scalar. Its explicit value will cancel. The expansion follows by Taylor expansion in normal coordinates on a fixed local neighbourhood, scaling coordinates by sqrt(epsilon), and bounding the complementary integral exponentially. Compact embeddedness ensures a uniform positive ambient separation for pairs outside a fixed intrinsic neighbourhood of the diagonal. Only integer powers survive Gaussian parity.

Construct

    q_N = rho^(-1/2) (1 + sum_{j=1}^N epsilon^j b_j).

The coefficient of b_j in the epsilon^j term of q_N T_epsilon(q_N rho) is 2b_j. All other terms involve already chosen b_1,...,b_{j-1}. Hence one recursively chooses smooth b_j such that

    q_N T_epsilon(q_N rho)=1+O(epsilon^(N+1))

in every fixed C^k norm. For sufficiently small epsilon, q_N is uniformly bounded above and below by positive constants. The first coefficient is

    b_1=-(L_1 sqrt(rho))/(2 sqrt(rho)).

Apply the stability lemma with exact v_epsilon, approximate q_N, and the positive-semidefinite kernel T_epsilon on L2(mu). (T_epsilon here acts on f rho when represented against mu.) Writing s=v_epsilon/q_N gives

    ||s-1||_L2 <= (sup s)^(1/2) O(epsilon^(N+1))
                 <= C epsilon^(N+1-d/8).

Use the direct kernel-smoothing upgrade in §3. More explicitly, for any prescribed finite expansion order J and derivative order k, take N sufficiently large; that estimate yields

    ||v_epsilon-q_N||_{C^k}=O(epsilon^(J+1)),

while q_N-q_J=O(epsilon^(J+1)) in C^k. Thus

    v_epsilon = rho^(-1/2)(1+epsilon b_1+...+epsilon^J b_J)
                  + O_{C^k}(epsilon^(J+1)).

This argument is not a circular use of a Laplace expansion on v_epsilon: Laplace expansion is applied only to q_N, whose coefficients are fixed smooth functions; a separate stability estimate and explicit kernel differentiation compare q_N with the exact solution.

## 5. Generator expansion

Substitute the scaling expansion into

    P_epsilon f=v_epsilon T_epsilon(v_epsilon rho f).

The first-order coefficient is

    (L_1(sqrt(rho) f))/sqrt(rho)
      - f (L_1 sqrt(rho))/sqrt(rho)
      = (1/4)(Delta_M f+grad(log rho).grad f).

The scalar E_M cancels exactly. Therefore

    P_epsilon=I-(epsilon/4)A+epsilon^2 B_2+...

as an asymptotic expansion on smooth functions, in every fixed output C^k norm. In particular,

    P_epsilon f=f-(epsilon/4)Af+O(epsilon^2).

## 6. Uniform gap and inverse bounds

The first scaling asymptotic yields constants 0<a<b such that a<=v_epsilon<=b. Hence

    <f,(I-P_epsilon)f>_mu
      = (1/2)(pi epsilon)^(-d/2)
        int int exp(-c(x,y)/epsilon) v_epsilon(x)v_epsilon(y)
                     (f(x)-f(y))^2 rho(x)rho(y) dxdy

dominates a constant times the unweighted local Gaussian form.

Let H_t be the ordinary heat semigroup on M with generator Delta_M. On a compact smooth closed manifold its heat kernel obeys a Gaussian upper bound

    h_t(x,y)<=C t^(-d/2) exp(-dist_M(x,y)^2/(C_0 t))

for small t. Set t=epsilon/C_0. Since c(x,y)=|x-y|^2<=dist_M(x,y)^2,

    (pi epsilon)^(-d/2) exp(-c(x,y)/epsilon)>=c_0 h_t(x,y).

Consequently the preceding form dominates c<f,(I-H_t)f>_{dx}. By the positive first nonzero Laplace eigenvalue of the connected compact manifold,

    <f,(I-H_t)f>_{dx} >= (1-exp(-lambda_1 t)) min_a int|f-a|^2 dx.

Uniform upper/lower density bounds compare the latter variance with Var_mu f. Since 1-exp(-lambda_1 t)>=c t for small t,

    <f,(I-P_epsilon)f>_mu >= c epsilon Var_mu f.

The Gaussian kernel is positive semidefinite; therefore 0<=P_epsilon<=I and, on mean-zero functions,

    ||P_epsilon||<=1-c epsilon,
    ||(I-P_epsilon^2)^(-1)||<=C/epsilon,
    ||R_epsilon||<=C.

The constants can be chosen uniformly over compact smooth families of strictly positive rho. The argument does not infer an inverse limit merely from pointwise generator convergence.

## 7. Resolvent/metric convergence

For smooth mean-zero h, let phi=A^(-1)h, so phi is smooth by elliptic regularity. The generator expansion gives

    (I-P_epsilon^2)phi=(epsilon/2)Aphi+O(epsilon^2)
                      =(epsilon/2)h+O(epsilon^2),
    (epsilon/2)P_epsilon h=(epsilon/2)h+O(epsilon^2).

Subtracting and applying the gap estimate gives

    ||R_epsilon h-phi||_L2
       <= (C/epsilon) O(epsilon^2)=O(epsilon).

Pair with h for the metric convergence. Since smooth mean-zero h are dense in L2_0(mu), and both R_epsilon and A^(-1) are uniformly bounded there, approximation proves strong convergence for every h in L2_0(mu).

## 8. Parameter derivatives and complete smooth resolvent expansion

For a smooth parameter family rho_theta, all previous constants and approximate-scaling coefficients are uniform on compact parameter sets. It remains to justify parameter derivatives, since uniform convergence alone does not imply their convergence.

At fixed epsilon, smooth dependence of the scaling on parameters follows from the implicit function theorem. In logarithmic variables, the derivative with respect to the scaling is I+P_epsilon. It is invertible on L2(mu) because P_epsilon is positive semidefinite. If a smooth right-hand side b is given, the L2 solution of (I+P_epsilon)w=b satisfies w=b-P_epsilon w, and the Gaussian kernel is smoothing, so w lies in every requested C^k space. This proves the required invertibility there and justifies differentiating the scaling equation.

Let w_i=partial_i log v_epsilon and h_i=partial_i log rho. Differentiating gives

    (I+P_epsilon) w_i=-P_epsilon h_i.

Since P_epsilon is positive semidefinite, (I+P_epsilon)^(-1) has L2(mu) norm at most one. Here is a direct comparison with the approximate scaling, which avoids assuming that parameter derivatives of the asymptotic expansion already hold.

Write q=q_N, r=q T_epsilon(q rho), and define the row-normalised approximate kernel

    P_q f = q T_epsilon(q rho f)/r.

Let w_q=partial_i log q and h=partial_i log rho. Differentiating the approximate row identity gives

    partial_i log r=(I+P_q)w_q+P_q h.

Subtracting from the exact identity yields

    (I+P_epsilon)(w_i-w_q)
      =-partial_i log r+(P_q-P_epsilon)(w_q+h).

The spatial scaling estimate already proves v_epsilon/q-1=O(epsilon^L) in any prescribed spatial C^k norm, for arbitrarily large L if N is chosen sufficiently large. Both kernels have positive bounded rows; their relative difference is of this order. The Schur test therefore bounds P_q-P_epsilon by O(epsilon^L) on L2(mu). Kernel differentiation gives the corresponding spatial C^k bounds with only polynomial losses in epsilon. Since partial_i log r=O(epsilon^(N+1)), the right-hand side is arbitrarily small in L2 and C^k after increasing N. Invert I+P_epsilon in L2, then use w_i-w_q=right-hand-side-P_epsilon(w_i-w_q) to gain the spatial C^k estimate by smoothing.

Repeated parameter differentiation gives the same leading operator I+P_epsilon multiplying the highest unknown logarithmic derivative. Every other term contains only already controlled lower derivatives, derivatives of the approximate scaling, or derivatives of its row residual. Induction proves the full scaling expansion in every mixed spatial/parameter C^k norm. This also yields the corresponding operator-kernel expansion and bounds for P_epsilon. Uniformity holds on compact parameter sets contained in the parameter domain.

To avoid a moving mean-zero subspace in the inverse problem, define

    Pi_theta f=(int f dmu_theta) 1,
    D_epsilon,theta=I-P_epsilon,theta^2+Pi_theta.

This operator acts as the identity on constants and as I-P_epsilon^2 on mean-zero functions. The gap estimate implies ||D_epsilon,theta^(-1)||<=C/epsilon uniformly, in L2(mu_theta) and hence in a fixed reference L2 space by uniform density equivalence. Parameter derivatives of P_epsilon and Pi_theta have bounded or at worst polynomial operator norms. The inverse identity gives, for each fixed multiindex alpha, a polynomial bound on partial_theta^alpha D_epsilon,theta^(-1); if derivatives of D are bounded, the simple bound C_alpha epsilon^(-|alpha|-1) suffices.

One can now construct a formal approximate solution

    phi_epsilon,N=phi_0+epsilon phi_1+...+epsilon^N phi_N

to the augmented equation

    D_epsilon,theta phi=(epsilon/2)P_epsilon h,

with residual O(epsilon^(N+2)) in every mixed C^k norm. At each order the unknown phi_j is acted on by A/2, so the recursion is an elliptic Poisson equation. Compatibility follows because every coefficient of the left-hand side and right-hand side has zero mu_theta mean; choose each phi_j to have zero mu_theta mean. The phi_j depend smoothly on theta by elliptic parameter theory.

Each coefficient is chosen with zero mu_theta mean. The exact solution has the same zero mean, since the right-hand side has zero mean and the augmented equation acts as identity on constants. After applying the exact inverse, the error e_N is O(epsilon^(N+1)) in L2; differentiating a fixed number of times costs only a polynomial power of epsilon. Taking N sufficiently large proves the expansion through any prescribed order in any prescribed parameter C^k norm.

The error also solves e_N=P_epsilon^2 e_N+residual, because its mean is zero. Gaussian kernel smoothing turns the L2 error into a spatial C^k error with a polynomial epsilon loss. Differentiating this identity in parameters, and using the already established L2 bounds for all requested parameter derivatives, proves the expansion in every mixed spatial/parameter C^k norm. Once more, an arbitrarily high residual order offsets every fixed polynomial loss.

Consequently G_epsilon=G_0+O_{C^k}(epsilon), including the explicit first-order term from §9a with O_{C^k}(epsilon^2) remainder. If the parameter map is a smooth immersion, its limiting Gram matrix is positive definite: a zero tangent metric implies its signed-density derivative is zero, hence the parameter tangent is zero. On a compact parameter subset the least metric eigenvalue has a positive lower bound. The ordinary finite-dimensional formulas for inverse metric, connection and curvature then yield their stated convergence.

## 9. Circle sanity check and sharpness

On a circle of radius R with uniform probability, write angular coordinate theta. The fixed ambient cost is

    c(theta,phi)=2R^2(1-cos(theta-phi)).

P_epsilon is convolution with the von Mises kernel, with Fourier eigenvalues

    p_k=I_k(kappa)/I_0(kappa), kappa=2R^2/epsilon.

For k!=0, A has eigenvalue k^2/R^2 and R_epsilon has eigenvalue

    r_epsilon,k=epsilon p_k/[2(1-p_k^2)].

The large-kappa expansion

    I_k(kappa) ~ exp(kappa)/sqrt(2pi kappa)
      [1-(4k^2-1)/(8kappa)
       +(4k^2-1)(4k^2-9)/(128kappa^2)+...]

gives

    p_k=1-k^2/(2kappa)+k^2(k^2-2)/(8kappa^2)+O(kappa^-3),
    log p_k=-k^2/(2kappa)-k^2/(4kappa^2)+O(kappa^-3).

Because r_epsilon,k=epsilon/[4sinh(-log p_k)],

    r_epsilon,k=R^2/k^2-epsilon/(4k^2)+O(epsilon^2).

Thus the O(epsilon) metric convergence rate is sharp even for a uniform, analytic density on a one-dimensional closed manifold. This first-order bias disappears for a periodised heat kernel, which is a different cost; one must not transfer its O(epsilon^2) rate to the fixed ambient-cost setting.

## 9a. Explicit general extrinsic-curvature correction

This is a stronger theorem than convergence alone. Write II for the second fundamental form, H=tr II for the (unnormalised) mean-curvature vector, and S_alpha for the shape operators in a local orthonormal normal frame. Set

    Q=sum_alpha S_alpha^2,
    <Qw,w>=sum_i |II(w,e_i)|^2.

This is independent of the chosen tangent/normal frames and is nonnegative. Let phi=A^(-1)h. Then the candidate sharp expansion is

    g_epsilon_mu(h mu,h mu)
      = int |grad phi|^2 dmu
        -(epsilon/4) int <Q grad phi,grad phi> dmu
        +O(epsilon^2).

It has been derived independently by root and scout_transport. The calculation is as follows.

Write T_epsilon=I+epsilon D+epsilon^2 F+..., with D=Delta/4+E. In the volume-space conjugation, sigma=sqrt(rho),

    sigma P_epsilon sigma^(-1)=q_epsilon T_epsilon q_epsilon,
    q_epsilon=v_epsilon sigma=1+epsilon b_1+epsilon^2 b_2+... .

The first coefficient L=D+2b_1 is sigma(-A/4)sigma^(-1). Subtracting L^2/2 from the second coefficient gives

    R + (2b_2-b_1^2)I,  where R=F-D^2/2.

Since the conjugated operator preserves sigma, this correction is

    R-(R sigma/sigma)I.

The operator R is formally self-adjoint. Its fourth-order principal symbol cancels because F has principal part Delta^2/32. A real formally self-adjoint differential operator cannot have a nonzero odd highest-order symbol, so R has order at most two. Therefore R=div(T grad)+W for a symmetric tensor T and scalar W. Determine T by applying R at a point x to f(exp_x v)=1/2(w.v)^2. In normal coordinates,

    |x-exp_x(v)|^2=|v|^2-|II(v,v)|^2/12+O(|v|^5),
    dvol(exp_x v)=[1-Ric(v,v)/6+O(|v|^3)]dv,
    Delta^2 f(x)=-(2/3)Ric(w,w).

Let Z be the Euclidean Gaussian of covariance I/2. The needed sixth moment is

    E[(w.Z)^2 |II(Z,Z)|^2]
      =(1/8)[|w|^2(|H|^2+2|II|^2)
             +4<H,II(w,w)>+8 sum_alpha |S_alpha w|^2].

Also E[(w.Z)^2 Ric(Z,Z)]=(1/4)[|w|^2 Scal+2Ric(w,w)]. Hence

    Ff(x)=E(x)|w|^2/4
           +<H,II(w,w)>/48
           +sum_alpha|S_alpha w|^2/24
           -Ric(w,w)/24,
    (D^2/2)f(x)=E(x)|w|^2/4-Ric(w,w)/48.

Subtract and use the contracted Gauss identity

    Ric(w,w)=<H,II(w,w)>-sum_alpha|S_alpha w|^2.

This yields

    <Tw,w>=Rf(x)=sum_alpha|S_alpha w|^2/16,

so T=Q/16. The scalar W cancels in R-(R sigma/sigma)I. Returning to L2(mu),

    P_epsilon=I-epsilon A/4
              +epsilon^2[A^2/32+C]+O(epsilon^3),
    C=(1/(16rho))div(rho Q grad).

The resolvent identity

    R_epsilon=(epsilon/4)[(I-P_epsilon)^(-1)-(I+P_epsilon)^(-1)]

or an approximate-Poisson-solution calculation now gives

    R_epsilon=A^(-1)+4epsilon A^(-1) C A^(-1)+O(epsilon^2)

on smooth mean-zero functions. Integration by parts proves the displayed metric correction. For a sphere of radius R, Q=R^(-2)I in every dimension, so the correction is -(epsilon/(4R^2))g_0, agreeing with the large-concentration von Mises–Fisher/Bessel calculation. For a circle this reduces to the explicit result of §9.

Potential significance: the leading bias depends on the embedding through its second fundamental form, even though the leading Wasserstein tangent metric depends only on the intrinsic geometry. It is nonpositive and yields an explicit first-order bias correction; this distinguishes the result sharply from heat-kernel Sinkhorn convergence, whose smooth-direction first-order term vanishes.

### A positive, computable second-order correction

The full smooth resolvent expansion also provides a correction that uses the entropic solution itself. Put phi_epsilon=R_epsilon h and define

    g_epsilon_corrected(h mu,h mu)
      =g_epsilon(h mu,h mu)
       +(epsilon/4)int <Q grad phi_epsilon,grad phi_epsilon> dmu.

Both summands are nonnegative, so this is a positive quadratic form. Since phi_epsilon=phi+O(epsilon) in C1 for smooth h,

    g_epsilon_corrected(h mu,h mu)=g_0(h mu,h mu)+O(epsilon^2).

For cross terms use polarisation. The result is uniform with derivatives over smooth finite-dimensional density families. This corrector requires the embedding curvature and the same entropic inverse used to compute the original metric. It does not require the unknown limiting Poisson solution phi. It has not been benchmarked against alternative algorithms or tested on real data.

## 10. Sources directly inspected

- Lavenant, H., Luckhardt, J., Mordant, G., Schmitzer, B., & Tamanini, L. The Riemannian geometry of Sinkhorn divergences. Annales de l'Institut Henri Poincare C 43 (2026), 1079–1144. https://doi.org/10.4171/AIHPC/165 ; https://arxiv.org/abs/2405.04987 . Full PDF read, especially §4.4, equation (35).
- Mordant, G., & Xu, L. (2026). Convergence of the Sinkhorn Riemannian metric for finitely supported measures. https://arxiv.org/abs/2609.14651 . Full PDF read, especially introduction and §3.3.
- Wormell, C. L., & Reich, S. (2021). Spectral convergence of diffusion maps: Improved error bounds and an alternative normalization. SIAM Journal on Numerical Analysis. https://doi.org/10.1137/20M1344093 ; https://arxiv.org/abs/2006.02037 . Full PDF read, especially Theorem 3.3, Theorem 3.5, Proposition 3.6. These give a related periodised-heat-kernel torus result, not the present fixed-cost proof.
- Marshall, N. F., & Coifman, R. R. (2019). Manifold learning with bi-stochastic kernels. IMA Journal of Applied Mathematics 84(3), 455–482. https://academic.oup.com/imamat/article-abstract/84/3/455/5289465 . See `novelty/source_inspection.md` for the reading record.

## 11. Boundaries

- No bounded-domain/reflection boundary layer is handled here.
- No unbounded Euclidean density, vanishing density, nonsmooth support, or singular base law is covered.
- Strong operator convergence and smooth finite-parameter geometry do not by themselves establish convergence of the global geodesic distances on unrestricted probability space.
- The small-bandwidth graph Laplacian limit is prior art. The candidate contribution is the rigorous inverse/Hessian limit, complete scaling error control, and continuous-density curvature convergence, plus sharp fixed-cost rate.
