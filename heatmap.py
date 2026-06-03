#!/usr/bin/env python3

import argparse
from pathlib import Path


REQUIRED_VARIABLES = ("x", "y", "z", "theta", "psi")
VARIABLE_CHOICES = REQUIRED_VARIABLES + ("psi_n",)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Sample ADIOS punctures onto a uniform grid and plot a heatmap."
    )
    parser.add_argument("adios_file", help="Input ADIOS BP file or BP directory")
    parser.add_argument(
        "-o",
        "--output",
        default="puncture_heatmap.png",
        help="Output PNG file (default: puncture_heatmap.png)",
    )
    parser.add_argument(
        "--x-var",
        choices=VARIABLE_CHOICES,
        default="y",
        help="ADIOS variable for the horizontal axis (default: y)",
    )
    parser.add_argument(
        "--y-var",
        choices=VARIABLE_CHOICES,
        default="z",
        help="ADIOS variable for the vertical axis (default: z)",
    )
    parser.add_argument(
        "--bins",
        nargs=2,
        type=int,
        metavar=("NX", "NY"),
        default=(256, 256),
        help="Number of uniform grid cells in x and y (default: 256 256)",
    )
    parser.add_argument(
        "--range",
        dest="histogram_range",
        nargs=4,
        type=float,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        help="Histogram range. Defaults to the finite data range.",
    )
    parser.add_argument(
        "--log",
        action="store_true",
        help="Use logarithmic color scaling for nonzero counts.",
    )
    parser.add_argument(
        "--equal-aspect",
        action="store_true",
        help="Use equal plot aspect ratio.",
    )
    parser.add_argument(
        "--cmap",
        default="viridis",
        help="Matplotlib colormap name (default: viridis)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the plot interactively after writing the PNG.",
    )
    return parser.parse_args()


def require_module(module_name, package_hint):
    try:
        return __import__(module_name)
    except ImportError as error:
        raise SystemExit(
            f"Missing Python module '{module_name}'. Install or load {package_hint} "
            "before running this script."
        ) from error


def read_required(stream, name, numpy):
    try:
        values = stream.read(name)
    except Exception as error:
        raise RuntimeError(f"Failed to read ADIOS variable '{name}'") from error
    return numpy.asarray(values).reshape(-1)


def variable_available(stream, name):
    available_variables = getattr(stream, "available_variables", None)
    if available_variables is None:
        return None

    try:
        return name in available_variables()
    except Exception:
        return None


def read_optional(stream, name, numpy):
    available = variable_available(stream, name)
    if available is False:
        return None

    try:
        values = stream.read(name)
    except Exception as error:
        if available is True:
            raise RuntimeError(f"Failed to read ADIOS variable '{name}'") from error
        return None

    return numpy.asarray(values).reshape(-1)


def append_step_arrays(stream, arrays, optional_names, optional_seen, numpy):
    for name in REQUIRED_VARIABLES:
        arrays[name].append(read_required(stream, name, numpy))

    for name in optional_names:
        values = read_optional(stream, name, numpy)
        present = values is not None
        previous = optional_seen.get(name)
        if previous is not None and previous != present:
            raise RuntimeError(f"ADIOS steps mix presence of optional variable '{name}'")
        optional_seen[name] = present
        if present:
            arrays.setdefault(name, []).append(values)


def finish_step_arrays(arrays, step_count, path, numpy):
    if step_count == 0:
        raise RuntimeError(f"No ADIOS steps found in {path}")

    return {
        name: numpy.concatenate(values)
        for name, values in arrays.items()
    }


def read_with_stream(adios2, path, optional_names, numpy):
    arrays = {name: [] for name in REQUIRED_VARIABLES}
    optional_seen = {}
    step_count = 0

    with adios2.Stream(str(path), "r") as stream:
        for _ in stream.steps():
            append_step_arrays(stream, arrays, optional_names, optional_seen, numpy)
            step_count += 1

    return finish_step_arrays(arrays, step_count, path, numpy)


def read_with_open(adios2, path, optional_names, numpy):
    arrays = {name: [] for name in REQUIRED_VARIABLES}
    optional_seen = {}
    step_count = 0

    with adios2.open(str(path), "r") as stream:
        for _ in stream:
            append_step_arrays(stream, arrays, optional_names, optional_seen, numpy)
            step_count += 1

    return finish_step_arrays(arrays, step_count, path, numpy)


def read_with_file_reader(adios2, path, optional_names, numpy):
    arrays = {}
    reader = adios2.FileReader(str(path))
    try:
        for name in REQUIRED_VARIABLES:
            arrays[name] = read_required(reader, name, numpy)

        for name in optional_names:
            values = read_optional(reader, name, numpy)
            if values is not None:
                arrays[name] = values
    finally:
        close = getattr(reader, "close", None)
        if close is not None:
            close()

    return arrays


def read_puncture_arrays(adios2, path, optional_names, numpy):
    if hasattr(adios2, "Stream"):
        return read_with_stream(adios2, path, optional_names, numpy)
    if hasattr(adios2, "open"):
        return read_with_open(adios2, path, optional_names, numpy)
    if hasattr(adios2, "FileReader"):
        return read_with_file_reader(adios2, path, optional_names, numpy)
    raise RuntimeError(
        "Unsupported ADIOS2 Python module: expected Stream, open, or FileReader"
    )


def selected_values(arrays, name):
    if name not in arrays:
        raise SystemExit(
            f"ADIOS variable '{name}' was requested but is not present in the file"
        )
    return arrays[name]


def finite_xy(x_values, y_values, numpy):
    mask = numpy.isfinite(x_values) & numpy.isfinite(y_values)
    return x_values[mask], y_values[mask], int(mask.size - numpy.count_nonzero(mask))


def histogram_range(args):
    if args.histogram_range is None:
        return None

    xmin, xmax, ymin, ymax = args.histogram_range
    if xmin >= xmax:
        raise SystemExit("--range requires XMIN < XMAX")
    if ymin >= ymax:
        raise SystemExit("--range requires YMIN < YMAX")

    return [[xmin, xmax], [ymin, ymax]]


def make_output_path(output):
    path = Path(output)
    parent = path.parent
    if parent != Path("."):
        parent.mkdir(parents=True, exist_ok=True)
    return path


def plot_heatmap(args, histogram, x_edges, y_edges, puncture_count, output_path, numpy):
    import matplotlib.pyplot as plt

    norm = None
    values = histogram.T
    colorbar_label = "punctures per cell"
    if args.log:
        from matplotlib.colors import LogNorm

        if histogram.max() > 0:
            values = numpy.ma.masked_where(values <= 0, values)
            norm = LogNorm(vmin=1, vmax=histogram.max())
            colorbar_label = "punctures per cell (log scale)"

    figure, axis = plt.subplots(figsize=(8.0, 6.5), constrained_layout=True)
    image = axis.imshow(
        values,
        origin="lower",
        extent=(x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]),
        aspect="auto",
        interpolation="nearest",
        cmap=args.cmap,
        norm=norm,
    )
    if args.equal_aspect:
        axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(args.x_var)
    axis.set_ylabel(args.y_var)
    axis.set_title(
        f"{args.y_var} vs {args.x_var} puncture density ({puncture_count} points)"
    )
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label(colorbar_label)
    figure.savefig(output_path, dpi=200)
    if args.show:
        plt.show()
    plt.close(figure)


def main():
    args = parse_args()
    if args.bins[0] <= 0 or args.bins[1] <= 0:
        raise SystemExit("--bins values must be positive")

    path = Path(args.adios_file)
    if not path.exists():
        raise SystemExit(f"ADIOS file does not exist: {path}")

    adios2 = require_module("adios2", "ADIOS2 Python bindings")
    numpy = require_module("numpy", "NumPy")
    matplotlib = require_module("matplotlib", "Matplotlib")
    if not args.show:
        matplotlib.use("Agg")

    optional_names = {
        name
        for name in (args.x_var, args.y_var)
        if name not in REQUIRED_VARIABLES
    }
    arrays = read_puncture_arrays(adios2, path, optional_names, numpy)
    x_values = selected_values(arrays, args.x_var)
    y_values = selected_values(arrays, args.y_var)
    if x_values.size != y_values.size:
        raise SystemExit(
            f"Selected variables have different lengths: "
            f"{args.x_var}={x_values.size}, {args.y_var}={y_values.size}"
        )

    x_values, y_values, dropped_count = finite_xy(x_values, y_values, numpy)
    if x_values.size == 0:
        raise SystemExit("No finite punctures remain after filtering")

    histogram, x_edges, y_edges = numpy.histogram2d(
        x_values,
        y_values,
        bins=args.bins,
        range=histogram_range(args),
    )
    output_path = make_output_path(args.output)
    plot_heatmap(args, histogram, x_edges, y_edges, x_values.size, output_path, numpy)

    print(f"Wrote {output_path}")
    print(f"Sampled {x_values.size} punctures into {args.bins[0]} x {args.bins[1]} bins")
    if dropped_count > 0:
        print(f"Dropped {dropped_count} non-finite punctures")


if __name__ == "__main__":
    main()
