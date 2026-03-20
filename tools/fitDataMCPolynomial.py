#Fit the data/mc polynomial for generating a warp and doing a warping study import argparse
import numpy as np
import uproot
from scipy.optimize import curve_fit

def polynomial(x, *coeffs):
    """Evaluate polynomial with given coefficients (p0 + p1*x + p2*x^2 + ...)"""
    return sum(c * x**i for i, c in enumerate(coeffs))

def get_histogram(filename, histname):
    """Extract bin centres, values and errors from a TH1D in a ROOT file."""
    with uproot.open(filename) as f:
        if histname not in f:
            raise KeyError(f"Histogram '{histname}' not found in {filename}. "
                           f"Available keys: {list(f.keys())}")
        hist = f[histname]
    
    values, edges = hist.to_numpy()
    centres = 0.5 * (edges[:-1] + edges[1:])
    errors = np.sqrt(hist.variances())
    
    return centres, values, errors

def compute_ratio(data_vals, data_errs, mc_vals):
    """
    Compute Data/MC ratio and propagate errors.
    Bins where MC is zero are masked out.
    """
    mask = mc_vals != 0
    ratio = np.where(mask, data_vals / mc_vals, np.nan)
    ratio_err = np.where(mask, data_errs / mc_vals, np.nan)
    return ratio, ratio_err, mask

def fit_polynomial(centres, ratio, ratio_err, degree, mask):
    """Fit a polynomial of given degree to the Data/MC ratio."""
    x = centres[mask]
    y = ratio[mask]
    yerr = ratio_err[mask]

    # Initial guess: all ones for constant-ish ratio near 1
    p0 = [1.0] + [0.0] * degree

    popt, pcov = curve_fit(
        polynomial, x, y,
        p0=p0,
        sigma=yerr,
        absolute_sigma=True
    )
    perr = np.sqrt(np.diag(pcov))
    return popt, perr

def print_results(var, degree, popt, perr):
    print(f"\nData/MC ratio polynomial fit for variable: {var}")
    print(f"Degree: {degree}")
    print("-" * 45)
    print(f"{'Parameter':<12} {'Value':>15} {'Error':>15}")
    print("-" * 45)
    for i, (val, err) in enumerate(zip(popt, perr)):
        print(f"  p{i:<9} {val:>15.6f} {err:>15.6f}")
    print("-" * 45)

    # Also print as a formatted expression
    terms = []
    for i, val in enumerate(popt):
        if i == 0:
            terms.append(f"{val:.6f}")
        elif i == 1:
            terms.append(f"({val:+.6f})*x")
        else:
            terms.append(f"({val:+.6f})*x^{i}")
    print(f"\nf(x) = {' + '.join(terms)}")

def main():
    parser = argparse.ArgumentParser(
        description="Fit a polynomial to the Data/MC ratio of a histogram."
    )
    parser.add_argument("var",      help="Variable name (histogram looked up as <var>_data / <var>_mc)")
    parser.add_argument("datafile", help="Path to ROOT file containing data histogram")
    parser.add_argument("mcfile",   help="Path to ROOT file containing MC histogram")
    parser.add_argument("degree",   type=int, help="Degree of polynomial to fit")
    args = parser.parse_args()

    hist_name_data = f"{args.var}_data"
    hist_name_mc   = f"{args.var}_mc"

    print(f"Loading '{hist_name_data}' from {args.datafile}")
    print(f"Loading '{hist_name_mc}'   from {args.mcfile}")

    centres, data_vals, data_errs = get_histogram(args.datafile, hist_name_data)
    _,       mc_vals,   _         = get_histogram(args.mcfile,   hist_name_mc)

    ratio, ratio_err, mask = compute_ratio(data_vals, data_errs, mc_vals)

    if mask.sum() < args.degree + 1:
        raise ValueError(
            f"Only {mask.sum()} valid bins but need at least {args.degree + 1} "
            f"to fit a degree-{args.degree} polynomial."
        )

    popt, perr = fit_polynomial(centres, ratio, ratio_err, args.degree, mask)
    print_results(args.var, args.degree, popt, perr)

if __name__ == "__main__":
    main()