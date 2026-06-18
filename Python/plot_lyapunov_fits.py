from __future__ import annotations

from pathlib import Path
import csv
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import Divider, Size


def build_sweep_arrays(treatment: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # FAK: 0 -> 0.4 -> 0.6 -> 1 (151 total)
    FAKs = np.concatenate((
        np.linspace(0, 0.4, 51),
        np.linspace(0.4, 0.6, 51)[1:],
        np.linspace(0.6, 1.0, 51)[1:],
    ))

    # AKT depends on treatment (151 total)
    if treatment == "DMSO":
        AKTs = np.concatenate((
            np.linspace(0, 0.6, 51),
            np.linspace(0.6, 1.0, 51)[1:],
            np.array([1.0] * 50),
        ))
    else:
        AKTs = np.concatenate((
            np.linspace(0, 0.3, 51),
            np.linspace(0.3, 0.6, 51)[1:],
            np.array([0.6] * 50),
        ))

    # ERK constant per treatment (repeated to 151)
    erk_val = 1.0 if treatment == "DMSO" else 0.2 if treatment == "ERL" else 0.1
    ERKs = np.array([erk_val] * 151, dtype=float)

    assert len(FAKs) == len(AKTs) == len(ERKs) == 151
    return AKTs, FAKs, ERKs


def fmt_param(v: float) -> str:
    v = float(np.round(v, 3))
    s = f"{v:.3f}".rstrip("0").rstrip(".")
    return s if s else "0"


def _summarize_indices_and_values(
    idx: np.ndarray,
    vals: np.ndarray,
    max_items: int = 20,
) -> tuple[str, str]:
    """Return (indices_str, values_str). If too many, truncates with a suffix."""
    if idx.size == 0:
        return "", ""

    take = min(max_items, idx.size)
    idx_take = idx[:take].astype(int)
    vals_take = vals[:take]

    idx_str = ";".join(map(str, idx_take.tolist()))
    val_str = ";".join([f"{float(v):.6g}" for v in vals_take.tolist()])

    if idx.size > max_items:
        idx_str += f";...(+{idx.size - max_items})"
        val_str += f";...(+{idx.size - max_items})"

    return idx_str, val_str


def mean_with_outlier_exclusion(
    arr: np.ndarray,
    *,
    method: str = "mad",   # "mad" or "iqr"
    mad_z: float = 3.5,    # common robust cutoff
    iqr_k: float = 1.5,    # Tukey fence multiplier
) -> dict:
    """
    Compute mean after excluding outliers, while tracking exclusions.
    Also excludes non-finite values (NaN/inf) and reports them.
    """
    a = np.asarray(arr).ravel()
    n_total = int(a.size)

    nonfinite_mask = ~np.isfinite(a)
    nonfinite_idx = np.where(nonfinite_mask)[0]
    nonfinite_vals = a[nonfinite_mask]

    finite_idx = np.where(~nonfinite_mask)[0]
    finite_vals = a[~nonfinite_mask]

    outlier_mask_finite = np.zeros(finite_vals.shape, dtype=bool)
    thresholds = ""

    if finite_vals.size > 0:
        if method.lower() == "mad":
            med = np.median(finite_vals)
            abs_dev = np.abs(finite_vals - med)
            mad = np.median(abs_dev)
            if mad > 0:
                robust_z = 0.6745 * (finite_vals - med) / mad
                outlier_mask_finite = np.abs(robust_z) > mad_z
                thresholds = f"mad_z>{mad_z}"
            else:
                thresholds = "mad=0 (no outliers)"
        elif method.lower() == "iqr":
            q1 = np.percentile(finite_vals, 25)
            q3 = np.percentile(finite_vals, 75)
            iqr = q3 - q1
            if iqr > 0:
                lo = q1 - iqr_k * iqr
                hi = q3 + iqr_k * iqr
                outlier_mask_finite = (finite_vals < lo) | (finite_vals > hi)
                thresholds = f"iqr_k={iqr_k} (lo={lo:.6g},hi={hi:.6g})"
            else:
                thresholds = "iqr=0 (no outliers)"
        else:
            raise ValueError(f"Unknown outlier method: {method!r}")

    outlier_idx = finite_idx[outlier_mask_finite]
    outlier_vals = a[outlier_idx]

    used_mask = np.ones_like(a, dtype=bool)
    used_mask[nonfinite_idx] = False
    used_mask[outlier_idx] = False
    used_vals = a[used_mask]

    n_nonfinite = int(nonfinite_idx.size)
    n_outlier = int(outlier_idx.size)
    n_used = int(used_vals.size)
    n_excluded_total = n_nonfinite + n_outlier

    mean_val = float(np.mean(used_vals)) if n_used > 0 else np.nan

    nonfinite_idx_str, nonfinite_vals_str = _summarize_indices_and_values(nonfinite_idx, nonfinite_vals)
    outlier_idx_str, outlier_vals_str = _summarize_indices_and_values(outlier_idx, outlier_vals)

    return {
        "mean": mean_val,
        "n_total": n_total,
        "n_used": n_used,
        "n_excluded_total": n_excluded_total,
        "n_nonfinite": n_nonfinite,
        "n_outlier": n_outlier,
        "nonfinite_indices": nonfinite_idx_str,
        "nonfinite_values": nonfinite_vals_str,
        "outlier_indices": outlier_idx_str,
        "outlier_values": outlier_vals_str,
        "outlier_method": method.lower(),
        "outlier_thresholds": thresholds,
    }


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def _interp_nans(y: np.ndarray) -> np.ndarray:
    """Linearly interpolate NaNs (and infs) so we can draw a continuous line."""
    y = np.asarray(y, dtype=float)
    x = np.arange(y.size)
    good = np.isfinite(y)
    if good.sum() == 0:
        return y.copy()
    if good.sum() == 1:
        out = np.full_like(y, np.nan)
        out[:] = y[good][0]
        return out
    return np.interp(x, x[good], y[good])


def smooth_series(y: np.ndarray, window: int = 9) -> np.ndarray:
    """
    Simple moving-average smoothing after NaN interpolation.
    - window should be odd and >= 3
    """
    y_filled = _interp_nans(y)

    w = int(window)
    if w < 3:
        return y_filled
    if w % 2 == 0:
        w += 1  # force odd

    pad = w // 2
    y_pad = np.pad(y_filled, (pad, pad), mode="edge")
    kernel = np.ones(w, dtype=float) / w
    y_smooth = np.convolve(y_pad, kernel, mode="valid")
    return y_smooth


def main() -> None:
    # Script is run from .../ECM_cell_cycle_model/Python
    script_dir = Path(__file__).resolve().parent
    lyap_root = script_dir / "LYAP_FIGURES"

    treatments = ["DMSO", "ERL", "LOW"]  # preserve order
    records: list[dict] = []

    # Outlier policy
    OUTLIER_METHOD = "mad"  # "mad" or "iqr"
    MAD_Z = 3.5
    IQR_K = 1.5

    for treatment in treatments:
        AKTs, FAKs, ERKs = build_sweep_arrays(treatment)

        for i in range(151):
            akt, fak, erk = AKTs[i], FAKs[i], ERKs[i]

            folder = f"AKT_{fmt_param(akt)}_FAK_{fmt_param(fak)}_ERK_{fmt_param(erk)}"
            run_dir = lyap_root / treatment / folder

            lyap_path = run_dir / "lyap.npy"
            ref_png = run_dir / "ref.png"

            # New marker files for classification (no JSON needed)
            steady_path = run_dir / "steady_lambda.npy"
            irreg_path = run_dir / "irregular.npy"

            # Defaults
            status = "ok"
            mean_val = np.nan
            dyn_state = "missing"
            steady_lambda_real = np.nan

            # Outlier/missing tracking fields (always present so CSV is consistent)
            n_total = n_used = n_excluded_total = 0
            n_nonfinite = n_outlier = 0
            nonfinite_indices = nonfinite_values = ""
            outlier_indices = outlier_values = ""
            outlier_method = OUTLIER_METHOD
            outlier_thresholds = ""

            # --- classify by file presence (preferred: lyap > steady > irregular) ---
            if lyap_path.exists():
                dyn_state = "limit_cycle"
            elif steady_path.exists():
                dyn_state = "steady_state"
                try:
                    lam = np.load(steady_path, allow_pickle=False)
                    # saved as complex; if NaN marker saved, this will become nan
                    steady_lambda_real = float(np.real(lam))
                except Exception:
                    steady_lambda_real = np.nan
            elif irreg_path.exists():
                dyn_state = "irregular_oscillation"
            else:
                dyn_state = "missing"

            # --- load lyap mean if present ---
            if lyap_path.exists():
                try:
                    arr = np.load(lyap_path, allow_pickle=False)

                    stats = mean_with_outlier_exclusion(
                        arr,
                        method=OUTLIER_METHOD,
                        mad_z=MAD_Z,
                        iqr_k=IQR_K,
                    )

                    mean_val = stats["mean"]
                    n_total = stats["n_total"]
                    n_used = stats["n_used"]
                    n_excluded_total = stats["n_excluded_total"]
                    n_nonfinite = stats["n_nonfinite"]
                    n_outlier = stats["n_outlier"]
                    nonfinite_indices = stats["nonfinite_indices"]
                    nonfinite_values = stats["nonfinite_values"]
                    outlier_indices = stats["outlier_indices"]
                    outlier_values = stats["outlier_values"]
                    outlier_method = stats["outlier_method"]
                    outlier_thresholds = stats["outlier_thresholds"]

                    if not np.isfinite(mean_val):
                        status = "invalid_mean (all excluded?)"
                    elif n_excluded_total > 0:
                        status = "ok_with_exclusions"
                    else:
                        status = "ok"

                except Exception as e:
                    status = f"load_error: {type(e).__name__}"
            else:
                # No lyap.npy: status based on classification markers / legacy artifacts
                if dyn_state == "steady_state":
                    status = "steady_state"
                elif dyn_state == "irregular_oscillation":
                    status = "irregular_oscillation"
                else:
                    if ref_png.exists():
                        status = "no_fit (ref.png exists, no marker files)"
                    else:
                        status = "missing (no lyap.npy, no ref.png, no markers)"

            records.append({
                "treatment": treatment,
                "index": i,
                "AKT": float(akt),
                "FAK": float(fak),
                "ERK": float(erk),
                "folder": folder,

                # primary outputs
                "mean_lyap": mean_val,
                "dyn_state": dyn_state,
                "steady_lambda_real": steady_lambda_real,
                "status": status,

                # paths
                "lyap_path": str(lyap_path),
                "steady_path": str(steady_path),
                "irreg_path": str(irreg_path),
                "run_dir": str(run_dir),

                # exclusion/missing diagnostics (only meaningful when lyap.npy exists)
                "n_total": n_total,
                "n_used": n_used,
                "n_excluded_total": n_excluded_total,
                "n_nonfinite": n_nonfinite,
                "n_outlier": n_outlier,
                "outlier_method": outlier_method,
                "outlier_thresholds": outlier_thresholds,
                "nonfinite_indices": nonfinite_indices,
                "nonfinite_values": nonfinite_values,
                "outlier_indices": outlier_indices,
                "outlier_values": outlier_values,
            })

    # Build plotting arrays in strict generation order (151 points per treatment)
    means_by_treat: dict[str, np.ndarray] = {}
    states_by_treat: dict[str, np.ndarray] = {}
    steady_by_treat: dict[str, np.ndarray] = {}

    for treatment in treatments:
        rows = [r for r in records if r["treatment"] == treatment]
        rows.sort(key=lambda r: r["index"])
        means_by_treat[treatment] = np.array([r["mean_lyap"] for r in rows], dtype=float)
        states_by_treat[treatment] = np.array([r["dyn_state"] for r in rows], dtype=object)
        steady_by_treat[treatment] = np.array([r["steady_lambda_real"] for r in rows], dtype=float)

    # CSV export: full results + issues-only subset
    out_dir = script_dir / "LYAP_FIGURES" / "_summaries"
    all_csv = out_dir / "lyap_means_with_markers.csv"
    issues_csv = out_dir / "lyap_issues_only.csv"

    fieldnames = [
        "treatment", "index", "AKT", "FAK", "ERK", "folder",
        "mean_lyap", "dyn_state", "steady_lambda_real", "status",
        "n_total", "n_used", "n_excluded_total", "n_nonfinite", "n_outlier",
        "outlier_method", "outlier_thresholds",
        "nonfinite_indices", "nonfinite_values",
        "outlier_indices", "outlier_values",
        "lyap_path", "steady_path", "irreg_path", "run_dir",
    ]
    write_csv(all_csv, records, fieldnames)

    issues = [
        r for r in records
        if (
            r["dyn_state"] != "limit_cycle"
            or (not np.isfinite(r["mean_lyap"]))
            or (r["n_excluded_total"] > 0)
            or ("load_error" in r["status"])
            or ("missing" in r["status"])
            or ("no_fit" in r["status"])
        )
    ]
    write_csv(issues_csv, issues, fieldnames)

    print(f"Wrote: {all_csv}")
    print(f"Wrote: {issues_csv}")
    print(f"Issues rows: {len(issues)} / {len(records)}")

    # Plot: treatment-colored Lyap lines; X markers for irregular; circles for steady eigenvalues
    params = {
        "axes.titley": 1,
        "axes.titlepad": -14,
        "axes.labelsize": 14,
        "axes.titlesize": 14,
        "axes.edgecolor": "black",
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    }
    plt.rcParams.update(params)

    # ===== SHARED axes-box geometry — keep IDENTICAL across matched figures =====
    BOX_W, BOX_H = 2.75, 2.75     # the plotting box, in inches (the thing being matched)
    LEFT, RIGHT  = 0.90, 0.15   # left margin holds the y-label + decimal y-ticks
    BOTTOM, TOP  = 0.55, 0.15
    figW = LEFT + BOX_W + RIGHT
    figH = BOTTOM + BOX_H + TOP

    fig = plt.figure(figsize=(figW, figH))
    fig.set_layout_engine("none")          # required — Divider places the axes manually
    h = [Size.Fixed(LEFT), Size.Fixed(BOX_W), Size.Fixed(RIGHT)]
    v = [Size.Fixed(BOTTOM), Size.Fixed(BOX_H), Size.Fixed(TOP)]
    div = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
    ax = fig.add_axes(div.get_position(), axes_locator=div.new_locator(nx=1, ny=1))

    color_map = {"DMSO": "C0", "ERL": "C1", "LOW": "C2"}
    x = np.linspace(0.0, 3.0, 151)

    SMOOTH_WINDOW = 11
    label_dict = {"ERL": "EGFRi", "LOW": "ERK = 0.1"}

    for treatment in treatments:
        y_raw = means_by_treat[treatment]
        states = states_by_treat[treatment]
        y_steady = steady_by_treat[treatment]
        c = color_map[treatment]

        # Smooth Lyap line (continuous even with NaNs by interpolating first)
        y_smooth = smooth_series(y_raw, window=SMOOTH_WINDOW)

        # --- NEW: do not draw line where there is a point attractor (steady_state) ---
        steady_mask = (states == "steady_state") & np.isfinite(y_steady)
        y_line = y_smooth.copy()
        y_line[steady_mask] = np.nan  # breaks the line at steady-state locations

        # Line only for non-steady regions
        ax.plot(
            x,
            y_line,
            linestyle="-",
            linewidth=3,
            color=c,
            label=label_dict.get(treatment, treatment),
        )

        # 'x' markers where oscillatory but irregular (no repeatable cycles)
        irreg = states == "irregular_oscillation"
        if np.any(irreg):
            ax.scatter(
                x[irreg],
                y_smooth[irreg],  # marker sits on the smoothed curve value
                marker="x",
                color=c,
                s=60,
                linewidths=1.8,
                alpha=0.9,
            )

        # Circles where steady-state eigenvalue exists (plot at eigenvalue, not on line)
        if np.any(steady_mask):
            ax.scatter(
                x[steady_mask],
                y_steady[steady_mask],
                marker="o",
                color=c,
                s=55,
                alpha=0.95,
            )

    # -- Annotation arrows pointing at notable features --
    # (x, y_tip, direction, treatment);  'up' = arrow points upward, 'down' = downward
    arrow_specs = [
        (0.05, -0.010, "up",   "LOW"),
        (0.75, -0.009, "up",   "ERL"),
        (1.00, -0.037, "up",   "DMSO"),
        (1.68, -0.018, "up",   "DMSO"),
        (1.98, -0.033, "up",   "ERL"),
        (2.30, -0.043, "up",   "LOW"),
        (2.45, -0.080, "down", "DMSO"),
        (2.98, -0.042, "down", "ERL"),
        (2.98, -0.055, "down", "DMSO"),
        (2.98, -0.060, "up",   "LOW"),
    ]
    ARROW_LEN = 0.012  # arrow length in y-units; tweak to taste

    for xa, y_tip, direction, treat in arrow_specs:
        dy = ARROW_LEN if direction == "up" else -ARROW_LEN
        ax.annotate(
            "",
            xy=(xa, y_tip),             # arrowhead sits at the feature
            xytext=(xa, y_tip - dy),    # tail offset
            arrowprops=dict(arrowstyle="-|>", color=color_map[treat],
                            lw=2, mutation_scale=11),
        )

    ax.set_xlim(0, 3)
    ax.set_xlabel("Stiffness Metric")
    ax.set_ylabel("Attractor Strength")

    # Get existing legend entries
    handles, labels = ax.get_legend_handles_labels()

    # Create new legend entries
    line_entry = Line2D([0], [0], color='black', lw=2)  # line
    circle_entry = Line2D([0], [0], marker='o', color='black', lw=0)  # circle
    x_entry = Line2D([0], [0], marker='x', color='black', lw=0, mew=2, ms=8)  # x

    # Append them
    handles.extend([line_entry, circle_entry, x_entry])
    labels.extend(['Cycle', 'Point', 'Indeterminate'])
    ax.legend(handles, labels, loc="lower left", fontsize=7, ncol=2,
              columnspacing=1, handlelength=1.4, labelspacing=0.3,
              borderpad=0.3, frameon=False)

    # Save as SVG (no bbox_inches='tight' -> the fixed box/figure size is preserved)
    out_svg = out_dir / "lyap_attractor_strength.svg"
    fig.savefig(out_svg, format="svg")
    print(f"Wrote: {out_svg}")

    plt.show()


if __name__ == "__main__":
    main()