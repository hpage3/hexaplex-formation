import argparse
import csv

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def d_to_q(d_angstrom):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 2.0 * np.pi / d_angstrom


def q_to_d(q_inv_angstrom):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 2.0 * np.pi / q_inv_angstrom


def two_theta_from_radius(radius_mm, detector_distance_mm):
    return np.arctan(radius_mm / detector_distance_mm)


def q_from_radius(radius_mm, detector_distance_mm, wavelength_angstrom):
    two_theta = two_theta_from_radius(radius_mm, detector_distance_mm)
    theta = two_theta / 2.0
    return 4.0 * np.pi * np.sin(theta) / wavelength_angstrom


def d_from_q(q_inv_angstrom):
    return q_to_d(q_inv_angstrom)


def radial_average(image, grid_limit, detector_distance_mm, wavelength_angstrom, bins):
    z_size, x_size = image.shape
    x_values = np.linspace(-grid_limit, grid_limit, x_size)
    z_values = np.linspace(-grid_limit, grid_limit, z_size)
    x_grid, z_grid = np.meshgrid(x_values, z_values)
    radius = np.sqrt(x_grid**2 + z_grid**2)

    r_min = 0.0
    r_max = float(radius.max())
    bin_edges = np.linspace(r_min, r_max, bins + 1)
    bin_indices = np.digitize(radius.ravel(), bin_edges, right=False) - 1
    bin_indices = np.where(bin_indices == bins, bins - 1, bin_indices)
    intensities = image.ravel()

    rows = []
    for index in range(bins):
        mask = bin_indices == index
        pixel_count = int(np.count_nonzero(mask))
        mean_intensity = float(np.mean(intensities[mask])) if pixel_count else np.nan
        r_min_mm = float(bin_edges[index])
        r_max_mm = float(bin_edges[index + 1])
        r_center_mm = (r_min_mm + r_max_mm) / 2.0
        two_theta_rad = float(
            two_theta_from_radius(r_center_mm, detector_distance_mm)
        )
        two_theta_deg = float(np.rad2deg(two_theta_rad))
        q_center = float(
            q_from_radius(r_center_mm, detector_distance_mm, wavelength_angstrom)
        )
        d_center = float(q_to_d(q_center)) if q_center != 0.0 else np.inf
        intensity_sum = float(np.sum(intensities[mask])) if pixel_count else np.nan
        rows.append(
            {
                "bin_index": index,
                "r_min_mm": r_min_mm,
                "r_max_mm": r_max_mm,
                "r_center_mm": r_center_mm,
                "radius_mm": r_center_mm,
                "two_theta_deg": two_theta_deg,
                "q_center_inv_angstrom": q_center,
                "q_Ainv": q_center,
                "d_center_angstrom": d_center,
                "d_A": d_center,
                "mean_intensity": mean_intensity,
                "intensity_mean": mean_intensity,
                "intensity_sum": intensity_sum,
                "intensity_total": intensity_sum,
                "pixel_count": pixel_count,
                "sample_count": pixel_count,
            }
        )

    return rows


def write_table(path, rows, delimiter=","):
    fieldnames = [
        "bin_index",
        "r_min_mm",
        "r_max_mm",
        "r_center_mm",
        "radius_mm",
        "two_theta_deg",
        "q_center_inv_angstrom",
        "q_Ainv",
        "d_center_angstrom",
        "d_A",
        "mean_intensity",
        "intensity_mean",
        "intensity_sum",
        "intensity_total",
        "pixel_count",
        "sample_count",
    ]
    with open(path, "w", encoding="utf-8", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def write_csv(path, rows):
    write_table(path, rows, delimiter=",")


def plot_profile(
    path,
    rows,
    plot_q_min=0.0,
    plot_log_y=False,
    normalize_plot=False,
    x_axis="q",
):
    plottable_rows = [
        row
        for row in rows
        if row["pixel_count"] > 0 and row["q_center_inv_angstrom"] >= plot_q_min
    ]
    if x_axis == "q":
        x_values = np.array([row["q_Ainv"] for row in plottable_rows])
        x_label = "q (Å^-1), convention q = 2π/d"
        title = "Radial diffraction profile in q-space"
    elif x_axis == "d":
        finite_rows = [row for row in plottable_rows if np.isfinite(row["d_A"])]
        x_values = np.array([row["d_A"] for row in finite_rows])
        plottable_rows = finite_rows
        x_label = "d-spacing (Å)"
        title = "Radial diffraction profile in d-spacing"
    elif x_axis == "radius":
        x_values = np.array([row["radius_mm"] for row in plottable_rows])
        x_label = "detector radius (mm)"
        title = "Radial diffraction profile by detector radius"
    else:
        raise ValueError(f"Unsupported x_axis: {x_axis}")

    mean_intensities = np.array([row["mean_intensity"] for row in plottable_rows])

    if normalize_plot and mean_intensities.size:
        max_intensity = np.max(mean_intensities)
        if max_intensity != 0:
            mean_intensities = mean_intensities / max_intensity

    if plot_log_y:
        positive = mean_intensities > 0
        x_values = x_values[positive]
        mean_intensities = mean_intensities[positive]

    plt.figure(figsize=(7, 4))
    plt.plot(x_values, mean_intensities, marker="o", markersize=3, linewidth=1)
    plt.xlabel(x_label)
    plt.ylabel("Normalized mean intensity" if normalize_plot else "Mean intensity")
    if plot_log_y:
        plt.yscale("log")
    if x_axis == "d":
        plt.gca().invert_xaxis()
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    return len(x_values)


def main():
    parser = argparse.ArgumentParser(
        description="Convert a 2D diffraction image into a 1D radial profile."
    )
    parser.add_argument("--input-npy", required=True, help="Input 2D .npy image.")
    parser.add_argument("--grid-limit", type=float, default=100.0)
    parser.add_argument("--detector-distance-mm", type=float, default=338.4)
    parser.add_argument("--wavelength-angstrom", type=float, default=0.7749)
    parser.add_argument("--bins", type=int, default=100)
    parser.add_argument("--output-prefix", default="outputs/radial_profile")
    parser.add_argument(
        "--table-format",
        choices=("csv", "tsv"),
        default="csv",
        help="Radial profile table format.",
    )
    parser.add_argument(
        "--include-empty-bins",
        action="store_true",
        help="Include empty radial bins in the CSV output; empty bins are never plotted.",
    )
    parser.add_argument(
        "--plot-log-y",
        action="store_true",
        help="Use a logarithmic y-axis in the PNG plot.",
    )
    parser.add_argument(
        "--plot-q-min",
        type=float,
        default=0.0,
        help="Minimum q value to include in the PNG plot only.",
    )
    parser.add_argument(
        "--normalize-plot",
        action="store_true",
        help="Normalize plotted intensities by the maximum plotted intensity.",
    )
    parser.add_argument(
        "--plot-space",
        choices=("q", "d", "radius", "both"),
        default="q",
        help="Plot q-space, d-spacing, detector radius, or both q and d plots.",
    )
    args = parser.parse_args()

    image = np.load(args.input_npy)
    if image.ndim != 2:
        raise ValueError(f"Expected a 2D input image, got shape {image.shape}")
    if args.bins < 1:
        raise ValueError("--bins must be positive")

    rows = radial_average(
        image,
        args.grid_limit,
        args.detector_distance_mm,
        args.wavelength_angstrom,
        args.bins,
    )
    non_empty_rows = [row for row in rows if row["pixel_count"] > 0]
    empty_bin_count = len(rows) - len(non_empty_rows)
    output_rows = rows if args.include_empty_bins else non_empty_rows

    table_extension = "tsv" if args.table_format == "tsv" else "csv"
    table_delimiter = "\t" if args.table_format == "tsv" else ","
    table_file = f"{args.output_prefix}.{table_extension}"
    write_table(table_file, output_rows, delimiter=table_delimiter)

    plot_outputs = []
    if args.plot_space == "both":
        for x_axis in ("q", "d"):
            png_file = f"{args.output_prefix}_{x_axis}.png"
            plotted_points = plot_profile(
                png_file,
                non_empty_rows,
                plot_q_min=args.plot_q_min,
                plot_log_y=args.plot_log_y,
                normalize_plot=args.normalize_plot,
                x_axis=x_axis,
            )
            plot_outputs.append((x_axis, png_file, plotted_points))
    else:
        png_file = f"{args.output_prefix}.png"
        plotted_points = plot_profile(
            png_file,
            non_empty_rows,
            plot_q_min=args.plot_q_min,
            plot_log_y=args.plot_log_y,
            normalize_plot=args.normalize_plot,
            x_axis=args.plot_space,
        )
        plot_outputs.append((args.plot_space, png_file, plotted_points))

    radial_means = np.array([row["mean_intensity"] for row in non_empty_rows])
    print(f"Input file: {args.input_npy}")
    print(f"Image shape: {image.shape}")
    print(f"Grid limit, mm: {args.grid_limit}")
    print(f"Number of bins: {args.bins}")
    print(f"Total radial bins requested: {len(rows)}")
    print(f"Non-empty radial bins: {len(non_empty_rows)}")
    empty_bin_action = "included" if args.include_empty_bins else "skipped"
    print(f"Empty radial bins {empty_bin_action}: {empty_bin_count}")
    print(f"Plot q minimum, inverse angstrom: {args.plot_q_min}")
    print(f"Plot log-y enabled: {args.plot_log_y}")
    print(f"Plot normalization enabled: {args.normalize_plot}")
    print(f"Plot space: {args.plot_space}")
    for x_axis, png_file, plotted_points in plot_outputs:
        print(f"Radial points plotted ({x_axis}): {plotted_points}")
        print(f"Output PNG ({x_axis}): {png_file}")
    print(f"Output table: {table_file}")
    print(f"2D input min intensity: {np.min(image)}")
    print(f"2D input max intensity: {np.max(image)}")
    print(f"2D input mean intensity: {np.mean(image)}")
    print(f"Radial mean min intensity: {np.min(radial_means)}")
    print(f"Radial mean max intensity: {np.max(radial_means)}")
    print(f"Radial mean mean intensity: {np.mean(radial_means)}")


if __name__ == "__main__":
    main()
