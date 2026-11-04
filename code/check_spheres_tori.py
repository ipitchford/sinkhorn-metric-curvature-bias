"""Assertion-free exact checks for the revised manuscript (B1, B2, B3 evidence).

1. Spheres S^m (m=1..4): Funk-Hecke eigenvalue ratio p_l = I_{l+nu}/I_nu checked against
   direct quadrature; first and second eps-coefficients of b_l checked at 60 digits.
2. Anisotropic flat torus S^1(r1) x S^1(r2) in R^4: exact factorised eigenvalues versus
   -(eps/4) <A^{-1} A_Q A^{-1}> with Q = diag(r1^-2, r2^-2).
3. Embedding dependence: round unit circle versus an ellipse of the same length 2*pi,
   uniform arclength density, phi = cos(s); predicted and observed first-order biases.
4. Circle operator-norm evidence: sup_k |b_k - R^2/k^2| / eps for decreasing eps.
Deterministic. Exits with status 1 if any check fails. Writes results/check_spheres_tori.json.
"""
import json, sys
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.special import ive, ellipe
from scipy.optimize import brentq

mp.mp.dps = 60
failures = []
out = {}

def check(name, cond, info):
    out.setdefault('checks', []).append({'name': name, 'pass': bool(cond), **info})
    if not cond:
        failures.append(name)

# ---------- 1. spheres ----------
def p_sphere(l, m, kappa):
    nu = mp.mpf(m - 1) / 2
    return mp.besseli(l + nu, kappa) / mp.besseli(nu, kappa)

def p_quad(l, m, kappa):
    # Funk-Hecke: ratio of int e^{kappa t} G_l(t)(1-t^2)^{(m-2)/2} dt, G_l normalised Gegenbauer
    if m == 1:
        f = lambda th: mp.exp(kappa * mp.cos(th))
        return mp.quad(lambda th: f(th) * mp.cos(l * th), [0, mp.pi]) / mp.quad(f, [0, mp.pi])
    lam = mp.mpf(m - 1) / 2
    def C(n, t):  # Gegenbauer recurrence, avoids hypergeometric evaluation
        c0, c1 = mp.mpf(1), 2 * lam * t
        if n == 0:
            return c0
        for j in range(1, n):
            c0, c1 = c1, (2 * t * (j + lam) * c1 - (j + 2 * lam - 1) * c0) / (j + 1)
        return c1
    G = lambda t: C(l, t) / C(l, mp.mpf(1))
    w = lambda t: mp.exp(kappa * t) * (1 - t * t)**(mp.mpf(m - 2) / 2)
    return mp.quad(lambda t: w(t) * G(t), [-1, 0, 1]) / mp.quad(w, [-1, 0, 1])

R = mp.mpf('1.3')
sph = []
for m in range(1, 5):
    for l in range(1, 4):
        q = p_quad(l, m, mp.mpf(5)); e = p_sphere(l, m, mp.mpf(5))
        check(f'funk_hecke_m{m}_l{l}', abs(q - e) < mp.mpf('1e-40'), {'quad_minus_bessel': mp.nstr(q - e, 5)})
        Lt = l * (l + m - 1); Lam = Lt / R**2
        c2_pred = ((m - 1)**2 + 2 * Lt - Lt**2 - 7) / (96 * R**4)
        ests = []
        for eps in [mp.mpf('1e-4'), mp.mpf('5e-5'), mp.mpf('2.5e-5')]:
            p = p_sphere(l, m, 2 * R**2 / eps)
            b = eps * p / (2 * (1 - p * p))
            c1 = (b * Lam - 1) / eps
            c2 = (b * Lam - 1 + eps / (4 * R**2)) / eps**2
            ests.append((c1, c2))
        # Richardson on c2 (error O(eps))
        c2_rich = 2 * ests[2][1] - ests[1][1]
        rec = {'m': m, 'l': l, 'c1_at_2.5e-5': mp.nstr(ests[2][0], 12), 'c1_pred': mp.nstr(-1 / (4 * R**2), 12),
               'c2_richardson': mp.nstr(c2_rich, 12), 'c2_pred': mp.nstr(c2_pred, 12)}
        sph.append(rec)
        check(f'sphere_c1_m{m}_l{l}', abs(ests[2][0] + 1 / (4 * R**2)) < 1e-3, rec)
        check(f'sphere_c2_m{m}_l{l}', abs(c2_rich - c2_pred) < 1e-5 * (1 + abs(c2_pred)), rec)
out['spheres'] = sph

# ---------- 2. flat torus in R^4 ----------
r1, r2 = mp.mpf(1), mp.mpf('0.6')
def p_circ(k, r, eps):
    kap = 2 * r**2 / eps
    return mp.besseli(k, kap) / mp.besseli(0, kap)
tor = []
for (k, l) in [(1, 0), (0, 1), (1, 1), (2, 1), (1, 2)]:
    Lam = k**2 / r1**2 + l**2 / r2**2
    LamQ = k**2 / r1**4 + l**2 / r2**4
    pred = -LamQ / (4 * Lam**2)
    vals = []
    for eps in [mp.mpf('1e-4'), mp.mpf('5e-5')]:
        p = p_circ(k, r1, eps) * p_circ(l, r2, eps)
        b = eps * p / (2 * (1 - p * p))
        vals.append((b - 1 / Lam) / eps)
    rich = 2 * vals[1] - vals[0]
    rec = {'k': k, 'l': l, 'slope_richardson': mp.nstr(rich, 12), 'predicted': mp.nstr(pred, 12)}
    tor.append(rec)
    check(f'flat_torus_{k}{l}', abs(rich - pred) < 1e-6, rec)
out['flat_torus_r1_1_r2_0.6'] = tor

# ---------- 3. circle vs ellipse of equal length ----------
a_ax = 1.3
b_ax = brentq(lambda b: 4 * a_ax * ellipe(1 - b * b / a_ax**2) - 2 * np.pi, 0.2, 1.3)
N = 1024
th = 2 * np.pi * np.arange(N) / N
X = np.stack([a_ax * np.cos(th), b_ax * np.sin(th)], 1)
J = np.sqrt((a_ax * np.sin(th))**2 + (b_ax * np.cos(th))**2)
kfreq = np.fft.fftfreq(N, 1 / N)
def d(f):
    return np.real(np.fft.ifft(1j * kfreq * np.fft.fft(f)))
Jh = np.fft.fft(J - J.mean()); Jh[0] = 0; Jh[1:] /= 1j * kfreq[1:]
s = th * J.mean() + np.real(np.fft.ifft(Jh))   # arclength, total length 2*pi
Lcheck = J.mean() * 2 * np.pi
phi = np.cos(s)
w = J * (2 * np.pi / N) / (2 * np.pi)          # uniform arclength probability
h = -d(d(phi) / J) / J
h -= np.dot(w, h)
kappa_c = a_ax * b_ax / J**3
phis = d(phi) / J
g0_ell = np.dot(w, phis**2); bias_ell = -0.25 * np.dot(w, kappa_c**2 * phis**2)
def sinkhorn_metric(eps):
    C = ((X[:, None, :] - X[None, :, :])**2).sum(-1)
    T = np.exp(-C / eps)
    v = np.ones(N) / np.sqrt(T @ w)
    for it in range(5000):
        r = v * (T @ (w * v)) - 1
        if np.max(np.abs(r)) < 1e-13:
            break
        v /= np.sqrt(1 + r)
    S = (np.sqrt(w) * v)[:, None] * T * (np.sqrt(w) * v)[None, :]
    lam, V = np.linalg.eigh(S)
    c = V.T @ (np.sqrt(w) * h)
    lam, c = lam[:-1], c[:-1]
    return float(eps / 2 * np.sum(c * c * lam / (1 - lam * lam))), float(np.max(np.abs(r)))
ell = []
for eps in [0.02, 0.01, 0.005]:
    g, res = sinkhorn_metric(eps)
    p = float(p_circ(1, mp.mpf(1), mp.mpf(eps)))
    g_circ = 0.5 * eps * p / (2 * (1 - p * p))  # h=cos s has |h|^2=1/2, eigen 1
    ell.append({'eps': eps, 'ellipse_slope': (g - g0_ell) / eps, 'circle_slope': (g_circ - 0.5) / eps,
                'scaling_residual': res})
rich_e = 2 * ell[2]['ellipse_slope'] - ell[1]['ellipse_slope']
rich_c = 2 * ell[2]['circle_slope'] - ell[1]['circle_slope']
out['equal_length_curves'] = {'ellipse_axes': [a_ax, b_ax], 'length_check': Lcheck, 'g0_ellipse': g0_ell,
                              'predicted_bias_ellipse': bias_ell, 'predicted_bias_circle': -0.125,
                              'richardson_ellipse': rich_e, 'richardson_circle': rich_c,
                              'shift_averaged_bias_ellipse': -0.125 * float(np.dot(w, kappa_c**2)), 'rows': ell}
check('equal_length_g0', abs(g0_ell - 0.5) < 1e-10, {'g0': g0_ell})
check('ellipse_bias', abs(rich_e - bias_ell) < 2e-3 * abs(bias_ell), {'rich': rich_e, 'pred': bias_ell})
check('circle_bias', abs(rich_c + 0.125) < 1e-4, {'rich': rich_c})
check('ellipse_differs_from_circle', abs(bias_ell + 0.125) > 1e-2, {'diff': bias_ell + 0.125})

# ---------- 4. circle operator norm ----------
rate = []
for eps in [1e-2, 3e-3, 1e-3, 3e-4, 1e-4]:
    k = np.arange(1, int(40 / np.sqrt(eps)) + 2)
    kap = 2 / eps
    p = ive(k, kap) / ive(0, kap)
    b = eps * p / (2 * (1 - p * p))
    bh = b + eps * k**2 / 4 * b * b
    dev = np.abs(b - 1.0 / k**2); devh = np.abs(bh - 1.0 / k**2)
    rate.append({'eps': eps, 'sup_over_eps': float(dev.max() / eps), 'argmax_k_sqrt_eps': float(k[dev.argmax()] * np.sqrt(eps)),
                 'corrected_sup_over_eps': float(devh.max() / eps), 'sup_over_sqrt_eps': float(dev.max() / np.sqrt(eps))})
out['circle_operator_norm'] = rate
check('circle_norm_bounded', max(r['sup_over_eps'] for r in rate) < 0.3, {})
out['failures'] = failures
path = Path(__file__).resolve().parents[1] / 'results' / 'check_spheres_tori.json'
path.write_text(json.dumps(out, indent=2, default=str) + '\n')
print(json.dumps({k: v for k, v in out.items() if k != 'checks'}, indent=1, default=str))
print('FAILURES:', failures)
sys.exit(1 if failures else 0)
