"""Compare Guo's Theorem 1 (Gaussian profile h(s)=exp(-s^2), bandwidth e^2=eps) with this paper.

Guo: (P-I)f/e^2 = C Lap f + E_in + E_ex + O(e^3) for the row-normalised P = Kf/K1, uniform density.
1. Symbolic moments for dimension d: on the unit S^d, Guo's eps^2 coefficient of p_l must equal
   Lam^2/32 - Lam/16, the coefficient from Prop 9.1 (Lam = l(l+d-1)); also compared with exact
   Bessel ratios at 50 digits for d=1..4.
2. Reduction: evaluating Guo's angular averages by the 15 pairings, for general M,
   second-order part = Lap^2/32 + (Q/16):Hess f + first-order terms (uses contracted Gauss).
   Checked symbolically in d, with the invariants as symbols.
3. Homogeneous example: product torus r=(1,0.6) in R^4, angular averages by quadrature, against
   the exact p_k p_l coefficient and against the Q form.
4. Non-homogeneous curve (ellipse): Guo's d=1 formula, including its gradient term 15/192 (k^2)' f',
   versus direct evaluation of the row-normalised operator.
5. First-order terms for every dimension (symbolic, using div(H.II)=gradS/2+divQ).
6. Torus of revolution: direct row-normalised coefficient versus the same prediction.
Exits 1 on any disagreement. Writes results/check_guo_comparison.json.
"""
import json, sys
from pathlib import Path
import mpmath as mp
import numpy as np
import sympy as sp

fails, out = [], {}
def check(name, ok, info=None):
    out[name] = {'pass': bool(ok), 'info': info}
    if not ok:
        fails.append(name)

d, Lam = sp.symbols('d Lambda', positive=True)
def mu0(k):                      # mu_k^(0) = Gamma((k+d)/2)/2, as ratio to mu_0^(0)
    return sp.rf(d / 2, sp.Rational(k, 2)) if k % 2 == 0 else None
def mu1(k):                      # mu_k^(1) = -(k+d-1) mu_{k-1}^(0), by parts (ratio to mu_0^(0))
    return -(k + d - 1) * mu0(k - 1)
m2, m4 = mu0(2), mu0(4)
C = sp.simplify(m2 / (2 * d))
check('C_is_quarter', sp.simplify(C - sp.Rational(1, 4)) == 0)

# 1. unit sphere S^d: Lap f=-Lam f, Ric=(d-1)g, S=d(d-1), |II(t,t)|^2=1
Lap, Lap2 = -Lam, Lam**2
trHessRic = (d - 1) * Lap
S = d * (d - 1)
Ein = m4 / (24 * d * (d + 2)) * (3 * Lap2 - 2 * trHessRic - 2 * Lap * S) + m2**2 / (12 * d**2) * Lap * S
Eex = -mu1(5) * (Lap / d) / 48 + mu1(3) * m2 * Lap / (48 * d)
guo_c2 = sp.simplify(Ein + Eex)
ours_c2 = Lam**2 / 32 - Lam / 16
check('sphere_symbolic_guo_equals_ours', sp.simplify(guo_c2 - ours_c2) == 0, str(sp.factor(guo_c2)))

mp.mp.dps = 50
rows = []
for dd in range(1, 5):
    nu = mp.mpf(dd - 1) / 2
    for l in (1, 2, 3):
        L = l * (l + dd - 1)
        vals = []
        for e in (mp.mpf('1e-4'), mp.mpf('5e-5')):
            kap = 2 / e
            p = mp.besseli(l + nu, kap) / mp.besseli(nu, kap)
            vals.append((p - 1 + L * e / 4) / e**2)
        rich = 2 * vals[1] - vals[0]
        pred = guo_c2.subs({d: dd, Lam: L})
        rows.append({'d': dd, 'l': l, 'exact': mp.nstr(rich, 10), 'guo': str(pred)})
        check(f'sphere_exact_d{dd}_l{l}', abs(rich - mp.mpf(sp.Rational(pred).p) / sp.Rational(pred).q) < 1e-5)
out['sphere_rows'] = rows

# 2. general reduction by pairings. Invariants at x, contracted with Hess f:
# avg |II(t,t)|^2 Hess(t,t) = [ I*Lf + 4*HII + 8*QH ]/(d(d+2)(d+4)), avg |II(t,t)|^2 = I/(d(d+2)),
# I=|H|^2+2|II|^2, HII=<H,II>:Hess, QH=Q:Hess.  Gauss: Ric:Hess = HII-QH, S=|H|^2-|II|^2.
I, Lf, HII, QH, Sc, GS, L2 = sp.symbols('I Lf HII QH Sc GS L2')
ric = HII - QH
Ein_g = m4 / (24 * d * (d + 2)) * (3 * L2 - 2 * ric - 2 * Lf * Sc - 3 * GS) + m2**2 / (12 * d**2) * Lf * Sc
Eex_g = -mu1(5) * (I * Lf + 4 * HII + 8 * QH) / (d * (d + 2) * (d + 4)) / 48 + mu1(3) * m2 * I / (d * (d + 2)) / (48 * d) * Lf
red = sp.simplify(Ein_g + Eex_g - (L2 / 32 + QH / 16 - GS / 32))
check('general_second_order_reduces_to_Q_over_16', red == 0, str(red))

# 3. product torus: Ric=0, S=0, d=2; f=cos(k a)cos(l b), arclength Hessian diag(-k^2/r1^2,-l^2/r2^2) f
r1, r2 = 1.0, 0.6
th = 2 * np.pi * (np.arange(4096) + .5) / 4096
c, s = np.cos(th), np.sin(th)
II2 = c**4 / r1**2 + s**4 / r2**2
tor = []
for k, l in [(1, 0), (0, 1), (1, 1), (2, 1)]:
    a, b = k**2 / r1**2, l**2 / r2**2
    hess = -a * c**2 - b * s**2
    lap = -(a + b)
    # Gaussian d=2: mu4/mu0=2, mu2/mu0=1, mu5^(1)/mu0=-6*2, mu3^(1)/mu0=-4*1, angular means
    guo = 2 / (24 * 2 * 4) * 3 * lap**2 + 12 / 48 * np.mean(II2 * hess) + (-4) * 1 * np.mean(II2) / (48 * 2) * (lap)
    qform = lap**2 / 32 + (-(a / r1**2 + b / r2**2)) / 16
    def co(n, r):
        return (n**2 / r**2) ** 2 / 32 - n**2 / r**4 / 16
    A, B = k**2 / r1**2 / 4, l**2 / r2**2 / 4
    exact = co(k, r1) + co(l, r2) + A * B
    tor.append({'k': k, 'l': l, 'guo': guo, 'Q_form': qform, 'exact': exact})
    check(f'torus_{k}{l}', abs(guo - exact) < 1e-10 and abs(qform - exact) < 1e-10, tor[-1])
out['product_torus'] = tor

# 4. ellipse: direct row-normalised operator at one point vs Guo d=1 and Q form
a_, b_ = 1.4, 0.8
N = 2048
t = 2 * np.pi * np.arange(N) / N
X, Y = a_ * np.cos(t), b_ * np.sin(t)
J = np.hypot(a_ * np.sin(t), b_ * np.cos(t))
freq = np.fft.fftfreq(N, 1 / N)
D = lambda g: np.real(np.fft.ifft(1j * freq * (np.abs(freq) <= 150) * np.fft.fft(g)))
ds = lambda g: D(g) / J
f = np.cos(t) + 0.5 * np.sin(2 * t)
k2 = (a_ * b_ / J**3)**2
f1 = ds(f); f2 = ds(f1); f4 = ds(ds(f2))
pred = f4 / 32 + k2 * f2 / 16 + (15 / 192) * ds(k2) * f1
ell = []
for i in (N // 8, N // 5, 3 * N // 7):
    vals = []
    for e in (0.01, 0.005):
        w = np.exp(-((X - X[i])**2 + (Y - Y[i])**2) / e) * J
        Pf = np.sum(w * f) / np.sum(w)
        vals.append((Pf - f[i] - e * f2[i] / 4) / e**2)
    rich = 2 * vals[1] - vals[0]
    ell.append({'i': i, 'richardson': rich, 'guo_Q_form': pred[i]})
    check(f'ellipse_{i}', abs(rich - pred[i]) < 2e-3 * (1 + abs(pred[i])), ell[-1])
out['ellipse'] = ell

# 5. first-order terms, every m: -gradS/32 + V/192 with V=gradI+4 div(H.II)+8 divQ and
# div(H.II)=gradS/2+divQ, against divQ/16 + gradE/4, E=I/48-S/12 (vector coefficients as symbols)
gI, gS, dQ = sp.symbols('gI gS dQ')
first = sp.simplify(-gS / 32 + (gI + 4 * (gS / 2 + dQ) + 8 * dQ) / 192 - (dQ / 16 + (gI / 48 - gS / 12) / 4))
check('first_order_every_m', first == 0, str(first))

# 6. torus of revolution (R=2,r=1): direct row-normalised eps^2 coefficient
u, v = sp.symbols('u v', real=True)
w = 2 + sp.cos(v)
gi = [1 / w**2, sp.Integer(1)]
k1, k2 = sp.cos(v) / w, sp.Integer(1)
ft = sp.cos(v) + sp.sin(u) * sp.sin(v) / 2 + sp.cos(2 * u) / 3
lap = lambda g: (sp.diff(w * gi[0] * sp.diff(g, u), u) + sp.diff(w * gi[1] * sp.diff(g, v), v)) / w
divQ = lambda g: (sp.diff(w * gi[0] * k1**2 * sp.diff(g, u), u) + sp.diff(w * k2**2 * sp.diff(g, v), v)) / w
Et = ((k1 + k2)**2 + 2 * (k1**2 + k2**2)) / 48 - 2 * k1 * k2 / 12
predt = lap(lap(ft)) / 32 + divQ(ft) / 16 + gi[1] * sp.diff(Et, v) * sp.diff(ft, v) / 4
Fp, F0, Fl = (sp.lambdify((u, v), e_, 'numpy') for e_ in (predt, ft, lap(ft)))
def X3(U, V):
    return (2 + np.cos(V)) * np.cos(U), (2 + np.cos(V)) * np.sin(U), np.sin(V)
def coeff(u0, v0, e, n=2400):
    g = 2 * np.pi * np.arange(n) / n
    UU, VV = np.meshgrid(g + u0, g + v0, indexing='ij')
    x, y, z = X3(UU, VV); a, b, c = X3(np.array(u0), np.array(v0))
    wt = np.exp(-((x - a)**2 + (y - b)**2 + (z - c)**2) / e) * (2 + np.cos(VV))
    return ((wt * F0(UU, VV)).sum() / wt.sum() - float(F0(u0, v0)) - e * float(Fl(u0, v0)) / 4) / e**2
trev = []
for u0, v0 in [(0.7, 0.9), (1.9, 2.2), (4.0, 0.4)]:
    c1, c2 = coeff(u0, v0, 0.01), coeff(u0, v0, 0.005)
    rich = 2 * c2 - c1
    rec = {'u': u0, 'v': v0, 'richardson': float(rich), 'predicted': float(Fp(u0, v0))}
    trev.append(rec)
    check(f'torus_rev_{u0}', abs(rich - rec['predicted']) < 2e-3 * (1 + abs(rec['predicted'])), rec)
out['torus_of_revolution_row'] = trev
out['failures'] = fails
Path(__file__).resolve().parents[1].joinpath('results', 'check_guo_comparison.json').write_text(json.dumps(out, indent=1, default=str) + '\n')
print(json.dumps({k: out[k] for k in ('torus_of_revolution_row', 'failures')}, indent=1, default=str))
for k in ('C_is_quarter','sphere_symbolic_guo_equals_ours','general_second_order_reduces_to_Q_over_16'): print(k, out[k])
sys.exit(1 if fails else 0)
