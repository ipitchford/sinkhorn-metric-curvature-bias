"""Deterministic check of the proposed fixed-cost curvature bias.

The output field corrected_error is the theoretical linear remainder
g_epsilon - g0 - epsilon*bias_coefficient. The separately named field
feasible_corrected_error uses the computable positive correction built from
the entropic potential. These are distinct quantities.
"""
import json
import argparse
from pathlib import Path
import platform
import numpy as np
import scipy
from scipy.linalg import solve


def run(n=768, epsilons=(0.12, 0.08, 0.04, 0.02, 0.01)):
    theta = 2 * np.pi * np.arange(n) / n
    rx, ry = 1.7, 0.8
    x = np.column_stack((rx*np.cos(theta), ry*np.sin(theta)))
    speed = np.sqrt(rx**2*np.sin(theta)**2 + ry**2*np.cos(theta)**2)
    dlogspeed = (rx**2-ry**2)*np.sin(theta)*np.cos(theta)/speed**2
    rho = np.exp(0.35*np.cos(theta) + 0.2*np.sin(2*theta))
    dlogrho = -0.35*np.sin(theta)+0.4*np.cos(2*theta)
    weights = rho*speed
    weights /= weights.sum()
    phi = np.sin(theta) + 0.3*np.cos(2*theta)
    dphi = np.cos(theta) - 0.6*np.sin(2*theta)
    ddphi = -np.sin(theta) - 1.2*np.cos(2*theta)
    score = -(ddphi+(dlogrho-dlogspeed)*dphi)/speed**2
    score -= weights@score
    curvature = rx*ry/speed**3
    g0 = float(weights@(dphi/speed)**2)
    bias_coefficient = float(-0.25*(weights@(curvature*dphi/speed)**2))
    distance2 = ((x[:,None,:]-x[None,:,:])**2).sum(axis=2)
    rootw = np.sqrt(weights)
    h = rootw*score
    projection = np.outer(rootw,rootw)
    identity = np.eye(n)
    results=[]
    for epsilon in epsilons:
        kernel=np.exp(-distance2/epsilon)
        u=np.ones(n)
        for step in range(1000):
            ku=kernel@(weights*u)
            residual=float(np.max(np.abs(u*ku-1)))
            if residual<2e-14:
                break
            u=np.sqrt(u/ku)
        if residual >= 2e-14:
            raise RuntimeError(
                f'Scaling did not converge for n={n}, epsilon={epsilon}: '
                f'residual={residual:.3e} after {step+1} iterations'
            )
        factor=rootw*u
        p=factor[:,None]*kernel*factor[None,:]
        left=solve(identity-p+projection,h,assume_a='pos')
        right=solve(identity+p,h,assume_a='pos')
        g=float(epsilon/4*(h@(left-right)))
        phie=epsilon/4*(left-right)/rootw
        frequencies=np.fft.fftfreq(n,1/n)
        dphie=np.fft.ifft(1j*frequencies*np.fft.fft(phie)).real
        feasible_corrected=float(g+epsilon/4*(weights@(curvature*dphie/speed)**2))
        err=g-g0-epsilon*bias_coefficient
        results.append(dict(epsilon=epsilon,metric=g,linear_bias=(g-g0)/epsilon,
                            corrected_error=err,corrected_error_over_epsilon2=err/epsilon**2,
                            feasible_corrected_metric=feasible_corrected,
                            feasible_corrected_error=feasible_corrected-g0,
                            feasible_corrected_error_over_epsilon2=(feasible_corrected-g0)/epsilon**2,
                            scaling_iterations=step,scaling_residual=residual))
    return dict(n=n,axes=[rx,ry],g0=g0,bias_coefficient=bias_coefficient,
                weighted_poisson_identity_error=float(weights@(score*phi)-g0),results=results)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parents[1]/'results')
    args=parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    all_results={'ellipse_bias_results':run()}
    for n in (768,1536):
        all_results[f'ellipse_bias_refinement_{n}']=run(n=n,epsilons=(.005,.0025))
    for name,answer in all_results.items():
        (args.output_dir/(name+'.json')).write_text(json.dumps(answer,indent=2)+'\n')
    metadata={
        'date':'2026-11-01','python':platform.python_version(),
        'numpy':np.__version__,'scipy':scipy.__version__,
        'method':'Deterministic periodic trapezoidal quadrature, symmetric Gaussian Sinkhorn scaling, dense positive-definite solves, FFT derivative for computable correction',
        'seed':'No randomness',
        'density':'rho(theta) proportional to exp(0.35*cos(theta)+0.2*sin(2*theta)) relative to arc length',
        'ellipse_axes':[1.7,0.8],
        'potential':'phi(theta)=sin(theta)+0.3*cos(2*theta)',
        'source_direction':'h=A phi, weighted mean removed to roundoff',
        'grid_sizes':[768,1536],
        'main_epsilons':[.12,.08,.04,.02,.01],
        'refinement_epsilons':[.005,.0025],
        'scaling_tolerance':2e-14,
        'threading':'Recommended OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1',
        'interpretation':'Numerical consistency check of an analytic theorem; no empirical data or stochastic inference.'
    }
    (args.output_dir/'ellipse_environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps({'g0':all_results['ellipse_bias_results']['g0'],
                      'predicted_first_order_coefficient':all_results['ellipse_bias_results']['bias_coefficient'],
                      'last_refinement':all_results['ellipse_bias_refinement_1536']['results'][-1],
                      'output_directory':str(args.output_dir)},indent=2))
