"""Symbolic certificate for the sphere expansion (Proposition on S^m).

Uses the large-argument expansion DLMF 10.40.1 of I_mu, truncated, to expand
  y_l = -log( I_{l+nu}(kappa) / I_nu(kappa) ),  nu=(m-1)/2,  kappa = 2R^2/eps,
and then b_l = eps/(4 sinh y_l).  Checks the closed forms stated in the manuscript:
  y_l = L u/2 + L u^2/4 - L(2(m-1)^2+4L-26) u^3/96 + O(u^4),  u=1/kappa,  L=l(l+m-1),
  b_l * L/R^2 = 1 - eps/(4R^2) + eps^2 ((m-1)^2+2L-L^2-7)/(96 R^4) + O(eps^3).
Exits with status 1 if any identity fails.  Writes results/sympy_sphere_coefficients.json.
"""
import json, sys
from pathlib import Path
import sympy as sp

u, mu, eps, R, m, l = sp.symbols('u mu epsilon R m l', positive=True)
ORDER = 4  # keep terms through u^ORDER

def trunc(expr, var, n):
    expr = sp.expand(expr)
    return sum(expr.coeff(var, j) * var**j for j in range(n + 1))

def a(k, mu_):
    p = sp.Integer(1)
    for j in range(1, k + 1):
        p *= (4 * mu_**2 - (2 * j - 1)**2)
    return p / (sp.factorial(k) * 8**k)

def log_series(mu_):
    X = sum((-1)**k * a(k, mu_) * u**k for k in range(1, ORDER + 1))
    out, P = 0, sp.Integer(1)
    for j in range(1, ORDER + 1):
        P = trunc(P * X, u, ORDER)
        out += sp.Rational((-1)**(j + 1), j) * P
    return trunc(out, u, ORDER)

nu = (m - 1) / 2
y = trunc(log_series(nu) - log_series(l + nu), u, ORDER)
L = l * (l + m - 1)
claims_y = {1: L / 2, 2: L / 4, 3: -L * (2 * (m - 1)**2 + 4 * L - 26) / 96}
ok = True
report = {}
for j, c in claims_y.items():
    d = sp.simplify(sp.expand(y.coeff(u, j) - c))
    report[f'y_coeff_u{j}_difference'] = str(d)
    ok &= (d == 0)

# b = eps/(4 sinh y); y in terms of eps up to eps^3 suffices for b through eps^2
ye = trunc(y.subs(u, eps / (2 * R**2)), eps, 3)
# 1/sinh(y) = (1/y)(1 - y^2/6 + 7 y^4/360 - ...); need relative accuracy eps^2
yt = sp.expand(ye / eps)          # y/eps, polynomial in eps
inv_lead = sp.series(1 / yt, eps, 0, 3).removeO()
corr = sp.series(1 - ye**2 / 6, eps, 0, 3).removeO()
b = trunc(sp.expand(inv_lead * corr / 4), eps, 2)
bL = sp.expand(b * L / R**2)
claims_b = {0: 1, 1: -1 / (4 * R**2), 2: ((m - 1)**2 + 2 * L - L**2 - 7) / (96 * R**4)}
for j, c in claims_b.items():
    d = sp.simplify(sp.expand(bL.coeff(eps, j) - c))
    report[f'bLambda_coeff_eps{j}_difference'] = str(d)
    ok &= (d == 0)
report['all_identities_hold'] = bool(ok)
out = Path(__file__).resolve().parents[1] / 'results' / 'sympy_sphere_coefficients.json'
out.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
sys.exit(0 if ok else 1)
