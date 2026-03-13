#!/usr/bin/env python3

import argparse
from contextlib import nullcontext
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pint import UnitRegistry
from scipy import constants as sp_constants

# -----------------------------
# Physical constants for Zimm plot
# -----------------------------
ureg = UnitRegistry()
N_SOLVENT = 1.331251
DN_DC = (0.37318 * ureg.milliliter / ureg.gram).to_base_units()
LAMBDA = (660 * ureg.nanometer).to_base_units()
AVOGADRO = sp_constants.Avogadro / ureg.mole
K_OPTICAL = ( 4.0 * np.pi**2 * N_SOLVENT**2 * DN_DC**2 ) / ( AVOGADRO * LAMBDA**4 )
K_OPTICAL = K_OPTICAL.to_base_units()
ANGLE_ROUND_DECIMALS = 1
CONC_COLOR_OVERRIDES = {
    0.0: "tab:blue",
    0.26893: "#f1c40f",
    0.79647: "tab:green",
    1.06413: "tab:red",
}

# -----------------------------
# Extract concentration from filename
# Expected format: protein_1p06413.dat
# -----------------------------
def extract_concentration(filename):
    # Look for pattern _NUMBERpNUMBER
    match = re.search(r'_(\d+p\d+)', filename)
    if not match:
        return None

    conc_str = match.group(1).replace('p', '.')
    return float(conc_str)  # in mg/mL


# -----------------------------
# Load two numeric columns
# -----------------------------

def load_two_columns(path):
    data = []
    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            try:
                angle = float(parts[0])
                intensity = float(parts[1])
                data.append((angle, intensity))
            except ValueError:
                continue
    if not data:
        raise ValueError(f"No numeric data found in {path}")
    arr = np.array(data)
    return arr[:, 0], arr[:, 1]


# -----------------------------
# Zimm plot helpers
# -----------------------------

def compute_baseline(df_mean):
    """Return baseline intensities (0 mg/mL) keyed by angle."""
    zero_mask = np.isclose(df_mean["concentration_mg_mL"], 0.0)
    baseline = df_mean[zero_mask].set_index("angle_deg")
    if baseline.empty:
        raise ValueError("Need a 0 mg/mL measurement to compute baselines.")
    return baseline


def compute_delta_df(df_mean, baseline):
    """Compute delta intensity relative to buffer for each concentration > 0."""
    positive = df_mean[df_mean["concentration_mg_mL"] > 0]
    rows = []

    for _, row in positive.iterrows():
        angle = row["angle_deg"]
        if angle not in baseline.index:
            continue
        buffer_intensity = baseline.at[angle, "intensity_mean"]
        buffer_sem = baseline.at[angle, "intensity_sem"]
        sample_sem = row["intensity_sem"]
        delta_value = row["intensity_mean"] - buffer_intensity
        delta_sem = np.sqrt(sample_sem**2 + buffer_sem**2)
        rows.append({
            "angle_deg": angle,
            "concentration_mg_mL": row["concentration_mg_mL"],
            "intensity_mean": row["intensity_mean"],
            "buffer_intensity": buffer_intensity,
            "intensity_sem": sample_sem,
            "buffer_sem": buffer_sem,
            "delta_intensity": delta_value,
            "delta_sem": delta_sem,
        })

    delta_df = pd.DataFrame(rows)
    if delta_df.empty:
        return delta_df
    return delta_df.sort_values(["concentration_mg_mL", "angle_deg"])


def compute_absolute_delta_df(df_mean):
    """Treat intensity itself as delta when no baseline is available."""
    positive = df_mean[df_mean["concentration_mg_mL"] > 0]
    rows = []
    for _, row in positive.iterrows():
        rows.append({
            "angle_deg": row["angle_deg"],
            "concentration_mg_mL": row["concentration_mg_mL"],
            "intensity_mean": row["intensity_mean"],
            "buffer_intensity": np.nan,
            "intensity_sem": row["intensity_sem"],
            "buffer_sem": np.nan,
            "delta_intensity": row["intensity_mean"],
            "delta_sem": row["intensity_sem"],
        })
    absolute_df = pd.DataFrame(rows)
    if absolute_df.empty:
        return absolute_df
    return absolute_df.sort_values(["concentration_mg_mL", "angle_deg"])


def compute_zimm_df(delta_df, zimm_constant):
    """Compute Zimm-plot-ready DataFrame for concentrations > 0."""
    if delta_df.empty:
        return pd.DataFrame()

    rows = []

    for _, row in delta_df.iterrows():
        angle = row["angle_deg"]
        delta_intensity = row["delta_intensity"]

        if np.isclose(delta_intensity, 0.0):
            continue

        conc_qty = (
            row["concentration_mg_mL"] * ureg.milligram / ureg.milliliter
        ).to_base_units()

        y_qty = (K_OPTICAL * conc_qty) / delta_intensity
        y_base = y_qty.to_base_units()
        y_sem_qty = (K_OPTICAL * conc_qty / (delta_intensity**2)) * row["delta_sem"]
        y_sem_base = y_sem_qty.to_base_units()

        inv_intensity = 1.0 / delta_intensity
        inv_intensity_sem = row["delta_sem"] / (delta_intensity**2)

        angle_rad = np.deg2rad(angle)
        sin_sq = np.sin(angle_rad / 2.0) ** 2
        x_value = sin_sq + zimm_constant * row["concentration_mg_mL"]

        rows.append({
            "angle_deg": angle,
            "concentration_mg_mL": row["concentration_mg_mL"],
            "delta_intensity": delta_intensity,
            "delta_sem": row["delta_sem"],
            "zimm_y_value": y_base.magnitude,
            "zimm_y_unit": str(y_base.units),
            "zimm_y_sem": y_sem_base.magnitude,
            "zimm_x_value": x_value,
            "sin_sq_half_angle": sin_sq,
            "inverse_intensity": inv_intensity,
            "inverse_intensity_sem": inv_intensity_sem,
        })

    zimm_df = pd.DataFrame(rows)
    if zimm_df.empty:
        return zimm_df

    return zimm_df.sort_values(["angle_deg", "zimm_x_value"])


def linear_regression(x_values, y_values):
    """Return slope and intercept for y = m*x + b; None if not enough variation."""
    x_arr = np.asarray(x_values, dtype=float)
    y_arr = np.asarray(y_values, dtype=float)
    if x_arr.size < 2:
        return None, None
    if np.allclose(x_arr, x_arr[0]):
        return None, None
    slope, intercept = np.polyfit(x_arr, y_arr, 1)
    return slope, intercept


def fit_zimm_plane(zimm_df):
    """Fit y = a + b*c + c*sin^2(theta/2) over all points."""
    if zimm_df.empty:
        return None
    conc = zimm_df["concentration_mg_mL"].to_numpy()
    sin_sq = zimm_df["sin_sq_half_angle"].to_numpy()
    y_vals = zimm_df["zimm_y_value"].to_numpy()
    ones = np.ones_like(conc)
    A = np.column_stack([ones, conc, sin_sq])
    coeffs, *_ = np.linalg.lstsq(A, y_vals, rcond=None)
    return {
        "intercept": coeffs[0],
        "slope_conc": coeffs[1],
        "slope_sin": coeffs[2],
        "max_conc": float(conc.max()),
        "max_sin": float(sin_sq.max()),
    }


def confirm_absolute_mode():
    """Ask user whether to continue without baseline data."""
    prompt = "No 0 mg/mL baseline found. Continue with absolute intensities only? [y/N]: "
    while True:
        try:
            reply = input(prompt).strip().lower()
        except EOFError:
            print("No input received; treating as 'no'.")
            return False
        if reply in ("y", "yes"):
            return True
        if reply in ("n", "no", ""):
            return False
        print("Please answer with 'y' or 'n'.")


# -----------------------------
# Main
# -----------------------------

def main():
    parser = argparse.ArgumentParser(description="Plot DLS Intensity vs Angle from folder.")
    parser.add_argument("folder", type=str, help="Folder containing .dat files")
    parser.add_argument("--out", type=str, default=None, help="Optional output PNG file")
    parser.add_argument(
        "--min-angle",
        type=float,
        default=None,
        help="Minimum scattering angle (deg) to include in tables and plots",
    )
    parser.add_argument(
        "--verbose", "---verbose",
        action="store_true",
        help="Print full (untruncated) DataFrames when reporting tables",
    )
    parser.add_argument(
        "--zimm-constant",
        type=float,
        default=100.0,
        help="Constant factor added in Zimm x-axis as sin^2(theta/2) + K*C (default: 100)",
    )
    args = parser.parse_args()

    folder = Path(args.folder)
    files = sorted(folder.glob("*.dat"))

    if not files:
        print("No .dat files found.")
        return

    all_data = []

    for file in files:
        conc = extract_concentration(file.name)
        
        if conc is None:
            print(f"Skipping {file.name} (no concentration found)")
            continue

        angles, intensities = load_two_columns(file)
        
        df_tmp = pd.DataFrame({
            "filename": file.name,
            "concentration_mg_mL": conc,
            "angle_deg": angles,
            "intensity": intensities
        })

        all_data.append(df_tmp)

    df = pd.concat(all_data, ignore_index=True)

    if args.min_angle is not None:
        df = df[df["angle_deg"] >= args.min_angle].copy()
        if df.empty:
            print(f"No measurements with angle >= {args.min_angle} deg. Nothing to plot.")
            return
        print(f"Filtering to scattering angles >= {args.min_angle} deg")

    df["angle_deg"] = df["angle_deg"].round(ANGLE_ROUND_DECIMALS)

    # Sort by concentration for readability
    df = df.sort_values(by="concentration_mg_mL")

    # Aggregate duplicate measurements per concentration and angle
    def _sem(series):
        if len(series) <= 1:
            return 0.0
        return series.std(ddof=1) / np.sqrt(len(series))

    df_mean = (
        df.groupby(["concentration_mg_mL", "angle_deg"], as_index=False)
        .agg(intensity_mean=("intensity", "mean"), intensity_sem=("intensity", _sem))
        .sort_values(["concentration_mg_mL", "angle_deg"])
    )

    zero_mask = np.isclose(df_mean["concentration_mg_mL"], 0.0)
    has_baseline = zero_mask.any()

    baseline = None
    delta_df = pd.DataFrame()
    zimm_df = pd.DataFrame()
    angle_regression_df = pd.DataFrame()
    plane_fit = None
    absolute_mode = False

    if has_baseline:
        baseline = compute_baseline(df_mean)
        delta_df = compute_delta_df(df_mean, baseline)
        zimm_df = compute_zimm_df(delta_df, args.zimm_constant)
        angle_rows = []
        for angle, group in zimm_df.groupby("angle_deg"):
            slope, intercept = linear_regression(group["concentration_mg_mL"], group["zimm_y_value"])
            if slope is None:
                continue
            angle_rows.append({
                "angle_deg": angle,
                "slope_dY_dConc": slope,
                "intercept_at_conc0": intercept,
            })
        if angle_rows:
            angle_regression_df = pd.DataFrame(angle_rows).sort_values("angle_deg")
        plane_fit = fit_zimm_plane(zimm_df)
    else:
        print("\nNo 0 mg/mL baseline detected.")
        if not confirm_absolute_mode():
            print("Aborting per user request.")
            return
        print("Continuing with absolute intensities only (baseline-dependent quantities are approximated).")
        absolute_mode = True
        delta_df = compute_absolute_delta_df(df_mean)
        if delta_df.empty:
            print("No positive concentration data available to build absolute Zimm plot.")
        else:
            zimm_df = compute_zimm_df(delta_df, args.zimm_constant)
            angle_rows = []
            for angle, group in zimm_df.groupby("angle_deg"):
                slope, intercept = linear_regression(group["concentration_mg_mL"], group["zimm_y_value"])
                if slope is None:
                    continue
                angle_rows.append({
                    "angle_deg": angle,
                    "slope_dY_dConc": slope,
                    "intercept_at_conc0": intercept,
                })
            if angle_rows:
                angle_regression_df = pd.DataFrame(angle_rows).sort_values("angle_deg")
            plane_fit = fit_zimm_plane(zimm_df)

    display_ctx = (
        pd.option_context(
            "display.max_rows", None,
            "display.max_columns", None,
            "display.width", None,
        )
        if args.verbose else
        nullcontext()
    )
    verbose_tag = "---verbose " if args.verbose else ""

    if args.verbose:
        with display_ctx:
            print(f"{verbose_tag}Aggregated intensity table:")
            print(df_mean)
            if not delta_df.empty:
                delta_title = "Delta intensity table (I_c - I_0)" if not absolute_mode else "Absolute intensity table (no baseline)"
                print(f"\n{verbose_tag}{delta_title}:")
                print(delta_df)
            elif has_baseline and not absolute_mode:
                print("\nDelta intensity table: no entries (missing buffer match).")

            if not zimm_df.empty:
                print(f"\n{verbose_tag}Zimm plot table:")
                print(zimm_df)
            elif has_baseline and not absolute_mode:
                print("\nZimm plot table: no entries (missing baseline match).")
            if not angle_regression_df.empty:
                print(f"\n{verbose_tag}Angle regression (y vs concentration) table:")
                print(angle_regression_df)
            if plane_fit:
                print(f"\n{verbose_tag}Global Zimm plane coefficients (y = a + b*c + c*sin^2):")
                plane_series = pd.Series({
                    "intercept": plane_fit["intercept"],
                    "slope_conc": plane_fit["slope_conc"],
                    "slope_sin": plane_fit["slope_sin"],
                })
                print(plane_series)
    else:
        print("Data tables suppressed; re-run with ---verbose to print full DataFrames.")
    # -----------------------------
    # Plot
    # -----------------------------

    unique_conc = sorted(df_mean["concentration_mg_mL"].unique())
    color_cycle = plt.rcParams.get("axes.prop_cycle", None)
    base_colors = []
    if color_cycle is not None:
        base_colors = color_cycle.by_key().get("color", [])
    if not base_colors:
        base_colors = [f"C{i}" for i in range(len(unique_conc))]
    conc_colors = {}
    auto_color_index = 0
    for conc in unique_conc:
        assigned = False
        for key, override in CONC_COLOR_OVERRIDES.items():
            if np.isclose(conc, key, atol=1e-4):
                conc_colors[float(conc)] = override
                assigned = True
                break
        if not assigned:
            conc_colors[float(conc)] = base_colors[auto_color_index % len(base_colors)]
            auto_color_index += 1

    fallback_color = base_colors[0] if base_colors else "black"

    def _color_for_conc(value):
        val = float(value)
        for key, override in CONC_COLOR_OVERRIDES.items():
            if np.isclose(val, key, atol=1e-4):
                return override
        return conc_colors.get(val, fallback_color)

    plt.figure(figsize=(8, 6))

    for conc, group in df_mean.groupby("concentration_mg_mL"):
        group = group.sort_values("angle_deg")
        plt.errorbar(
            group["angle_deg"],
            group["intensity_mean"],
            yerr=group["intensity_sem"],
            marker='o',
            linestyle='none',
            capsize=3,
            color=_color_for_conc(conc),
            label=f"{conc:.5f} mg/mL"
        )

    plt.xlabel("Scattering Angle (deg)")
    plt.ylabel("Intensity (a.u.)")
    plt.title("Intensity vs Angle")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    if args.out:
        plt.savefig(args.out, dpi=300)
        print(f"Plot saved to {args.out}")
        plt.close()
    else:
        plt.show()

    if not delta_df.empty:
        plt.figure(figsize=(8, 6))
        for conc, group in delta_df.groupby("concentration_mg_mL"):
            group = group.sort_values("angle_deg")
            plt.errorbar(
                group["angle_deg"],
                group["delta_intensity"],
                yerr=group["delta_sem"],
                marker='o',
                linestyle='none',
                capsize=3,
                color=_color_for_conc(conc),
                label=f"{conc:.5f} mg/mL"
            )

        if not absolute_mode:
            plt.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        plt.xlabel("Scattering Angle (deg)")
        delta_ylabel = "ΔI = I_c - I_0 (a.u.)" if not absolute_mode else "Intensity (a.u.)"
        delta_title = "ΔI vs Angle" if not absolute_mode else "Intensity vs Angle (no baseline)"
        plt.ylabel(delta_ylabel)
        plt.title(delta_title)
        plt.legend(title="Concentration (mg/mL)")
        plt.grid(True)
        plt.tight_layout()

        if args.out:
            delta_out_path = Path(args.out)
            delta_file = delta_out_path.with_name(f"{delta_out_path.stem}_delta{delta_out_path.suffix}")
            plt.savefig(delta_file, dpi=300)
            print(f"ΔI plot saved to {delta_file}")
            plt.close()
        else:
            plt.show()
    else:
        print("Skipping ΔI plot (no valid data).")

    if not zimm_df.empty:
        plt.figure(figsize=(8, 6))
        for conc, group in zimm_df.groupby("concentration_mg_mL"):
            group = group.sort_values("zimm_x_value")
            color = _color_for_conc(conc)
            plt.errorbar(
                group["zimm_x_value"],
                group["zimm_y_value"],
                yerr=group["zimm_y_sem"],
                marker='o',
                linestyle='none',
                capsize=3,
                color=color,
                label=f"{conc:.5f} mg/mL"
            )
            plt.plot(
                group["zimm_x_value"],
                group["zimm_y_value"],
                linestyle="-",
                alpha=0.6,
                color=color,
            )

        if plane_fit:
            sin_vals = np.linspace(0, plane_fit["max_sin"], 200)
            conc_vals = np.linspace(0, plane_fit["max_conc"], 200)
            y_conc_zero = plane_fit["intercept"] + plane_fit["slope_sin"] * sin_vals
            x_conc_zero = sin_vals  # concentration term drops out
            plt.plot(
                x_conc_zero,
                y_conc_zero,
                linestyle="--",
                linewidth=1.4,
                color="black",
                label="θ sweep (c→0)",
            )
            y_angle_zero = plane_fit["intercept"] + plane_fit["slope_conc"] * conc_vals
            x_angle_zero = args.zimm_constant * conc_vals
            plt.plot(
                x_angle_zero,
                y_angle_zero,
                linestyle=":",
                linewidth=1.4,
                color="black",
                label="c sweep (θ→0)",
            )

        y_unit = zimm_df["zimm_y_unit"].iloc[0]
        plt.xlabel(f"sin^2(theta/2) + {args.zimm_constant:g} * concentration (mg/mL)")
        plt.ylabel(f"K*c/(I - I0) [{y_unit}]")
        plt.title("Zimm Plot")
        plt.legend(title="Concentration (mg/mL)")
        plt.grid(True)
        plt.tight_layout()

        if args.out:
            zimm_out_path = Path(args.out)
            zimm_file = zimm_out_path.with_name(f"{zimm_out_path.stem}_zimm{zimm_out_path.suffix}")
            plt.savefig(zimm_file, dpi=300)
            print(f"Zimm plot saved to {zimm_file}")
            plt.close()
        else:
            plt.show()

        # Additional plot: inverse intensity vs x-axis (no K*c scaling)
        plt.figure(figsize=(8, 6))
        for conc, group in zimm_df.groupby("concentration_mg_mL"):
            group = group.sort_values("zimm_x_value")
            color = _color_for_conc(conc)
            plt.errorbar(
                group["zimm_x_value"],
                group["inverse_intensity"],
                yerr=group["inverse_intensity_sem"],
                marker='o',
                linestyle='none',
                capsize=3,
                color=color,
                label=f"{conc:.5f} mg/mL"
            )
            plt.plot(
                group["zimm_x_value"],
                group["inverse_intensity"],
                linestyle="-",
                alpha=0.6,
                color=color,
            )

        plt.xlabel(f"sin^2(theta/2) + {args.zimm_constant:g} * concentration (mg/mL)")
        plt.ylabel("1 / (I - I0) [1/a.u.]")
        plt.title("Inverse Intensity vs sin^2(theta/2) + K * concentration")
        plt.legend(title="Concentration (mg/mL)")
        plt.grid(True)
        plt.tight_layout()

        if args.out:
            zimm_inv_file = zimm_out_path.with_name(f"{zimm_out_path.stem}_zimm_inverse{zimm_out_path.suffix}")
            plt.savefig(zimm_inv_file, dpi=300)
            print(f"Inverse-intensity plot saved to {zimm_inv_file}")
            plt.close()
        else:
            plt.show()
    else:
        print("Skipping Zimm plot (no valid data).")

    print("\nAggregated DataFrame preview:")
    print(df_mean.head())

    if not delta_df.empty:
        print("\nΔI DataFrame preview:")
        print(delta_df.head())

    if not zimm_df.empty:
        print("\nZimm DataFrame preview:")
        print(zimm_df.head())

    if not args.verbose and not angle_regression_df.empty:
        print("\nAngle regression summary (y vs concentration):")
        print(angle_regression_df)

    if not args.verbose and plane_fit:
        print("\nGlobal Zimm plane coefficients (y = a + b*c + c*sin^2):")
        print(f"intercept={plane_fit['intercept']:.6e}, slope_conc={plane_fit['slope_conc']:.6e}, slope_sin={plane_fit['slope_sin']:.6e}")


if __name__ == "__main__":
    main()
