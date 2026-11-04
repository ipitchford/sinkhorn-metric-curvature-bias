"""Non-umbilic surface check: torus of revolution in R^3 at nonconstant density.

x(u,v) = ((a+b cos v) cos u, (a+b cos v) sin u, b sin v), a=2, b=1; density on area
proportional to exp(0.3 cos v + 0.2 sin v).  Rotational symmetry reduces the Sinkhorn
operator on u-Fourier mode n to a 1-D operator in v with exact kernel
  K_n(v,v') = ive(n, 2 r r'/eps) exp(-((r-r')^2 + b^2 (sin v - sin v')^2)/eps),  r=a+b cos v.
Directions phi = cos(n u) psi(v), n=0,1.  Principal curvatures cos v / r (parallel) and 1/b
(meridian), so Q is not a multiple of the identity and Ric = Gauss curvature varies in sign.
Compares (g_eps-g_0)/eps (Richardson) with -(1/4) int <Q grad phi, grad phi> dmu.
Deterministic; exits 1 on failure; writes results/check_torus_revolution.json.
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.special import ive

a, b = 2.0, 1.0
def run(N, eps, n, psi_fun):
    v = 2 * np.pi * np.arange(N) / N
    r = a + b * np.cos(v)
    rho = np.exp(0.3 * np.cos(v) + 0.2 * np.sin(v))
    om = rho * r; om /= om.sum()                     # v-marginal of mu
    kf = np.fft.fftfreq(N, 1 / N)
    d = lambda f: np.real(np.fft.ifft(1j * kf * np.fft.fft(f)))
    psi = psi_fun(v); dpsi = d(psi)
    H = n * n * psi / r**2 - d(rho * r * dpsi) / (rho * r * b * b)
    if n == 0:
        H = H - np.dot(om, H)
    g0 = np.dot(om, n * n * psi**2 / r**2 + dpsi**2 / b**2)
    kap_u = np.cos(v) / r; kap_v = 1 / b
    bias = -0.25 * np.dot(om, kap_u**2 * n * n * psi**2 / r**2 + kap_v**2 * dpsi**2 / b**2)
    R1, R2 = np.meshgrid(r, r, indexing='ij'); S1, S2 = np.meshgrid(np.sin(v), np.sin(v), indexing='ij')
    E = np.exp(-((R1 - R2)**2 + b * b * (S1 - S2)**2) / eps)
    K0 = ive(0, 2 * R1 * R2 / eps) * E
    s = 1 / np.sqrt(K0 @ om)
    for it in range(20000):
        res = s * (K0 @ (om * s)) - 1
        if np.max(np.abs(res)) < 1e-13:
            break
        s /= np.sqrt(1 + res)
    Kn = K0 if n == 0 else ive(n, 2 * R1 * R2 / eps) * E
    Sm = (np.sqrt(om) * s)[:, None] * Kn * (np.sqrt(om) * s)[None, :]
    lam, V = np.linalg.eigh(Sm)
    c = V.T @ (np.sqrt(om) * H)
    if n == 0:
        lam, c = lam[:-1], c[:-1]
    g = eps / 2 * np.sum(c * c * lam / (1 - lam * lam))
    return dict(g=float(g), g0=float(g0), bias=float(bias), residual=float(np.max(np.abs(res))), iters=it)

cases = {0: lambda v: np.sin(v) + 0.3 * np.cos(2 * v), 1: lambda v: 1 + 0.4 * np.sin(v)}
out, fails = {}, []
for n, f in cases.items():
    rows = []
    for eps in [0.02, 0.01, 0.005]:
        rr = run(512, eps, n, f); rr['eps'] = eps; rr['slope'] = (rr['g'] - rr['g0']) / eps; rows.append(rr)
    ref = run(768, 0.005, n, f)
    rich = 2 * rows[2]['slope'] - rows[1]['slope']
    rec = dict(rows=rows, richardson=rich, predicted=rows[0]['bias'], refinement_768_diff=abs(ref['g'] - rows[2]['g']))
    out[f'mode_{n}'] = rec
    if not (abs(rich - rec['predicted']) < 0.01 * abs(rec['predicted']) and rec['refinement_768_diff'] < 1e-9
            and max(x['residual'] for x in rows) < 1e-12):
        fails.append(n)
out['failures'] = fails
Path(__file__).resolve().parents[1].joinpath('results', 'check_torus_revolution.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=1)); sys.exit(1 if fails else 0)
