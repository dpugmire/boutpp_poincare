#!/usr/bin/env python3

import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot Poincare punctures from an ADIOS BP file."
    )
    parser.add_argument("adios_file", help="Input ADIOS BP file or BP directory")
    parser.add_argument(
        "output_name",
        help="Output basename. Writes <output_name>_yz.png and <output_name>_theta_psi.png",
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


def read_with_stream(adios2, path, numpy):
    arrays = {"y": [], "z": [], "theta": [], "psi": []}
    step_count = 0

    with adios2.Stream(str(path), "r") as stream:
        for _ in stream.steps():
            for name in arrays:
                arrays[name].append(read_required(stream, name, numpy))
            step_count += 1

    if step_count == 0:
        raise RuntimeError(f"No ADIOS steps found in {path}")
    return {name: numpy.concatenate(values) for name, values in arrays.items()}


def read_with_open(adios2, path, numpy):
    arrays = {"y": [], "z": [], "theta": [], "psi": []}
    step_count = 0

    with adios2.open(str(path), "r") as stream:
        for _ in stream:
            for name in arrays:
                arrays[name].append(read_required(stream, name, numpy))
            step_count += 1

    if step_count == 0:
        raise RuntimeError(f"No ADIOS steps found in {path}")
    return {name: numpy.concatenate(values) for name, values in arrays.items()}


def read_with_file_reader(adios2, path, numpy):
    arrays = {}
    reader = adios2.FileReader(str(path))
    try:
        for name in ("y", "z", "theta", "psi"):
            try:
                arrays[name] = numpy.asarray(reader.read(name)).reshape(-1)
            except Exception as error:
                raise RuntimeError(f"Failed to read ADIOS variable '{name}'") from error
    finally:
        close = getattr(reader, "close", None)
        if close is not None:
            close()

    return arrays


def read_puncture_arrays(adios2, path, numpy):
    if hasattr(adios2, "Stream"):
        return read_with_stream(adios2, path, numpy)
    if hasattr(adios2, "open"):
        return read_with_open(adios2, path, numpy)
    if hasattr(adios2, "FileReader"):
        return read_with_file_reader(adios2, path, numpy)
    raise RuntimeError(
        "Unsupported ADIOS2 Python module: expected Stream, open, or FileReader"
    )


def output_paths(output_name):
    base = Path(output_name)
    parent = base.parent
    if parent != Path("."):
        parent.mkdir(parents=True, exist_ok=True)

    stem = base.name
    if base.suffix:
        stem = base.stem

    yz = parent / f"{stem}_yz.png"
    theta_psi = parent / f"{stem}_theta_psi.png"
    return yz, theta_psi


def plot_scatter(matplotlib, x, y, xlabel, ylabel, title, output_path,
                 equal_aspect=False):
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(7.0, 6.0), constrained_layout=True)
    axis.scatter(x, y, s=2.0, alpha=0.65, linewidths=0)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.set_title(title)
    axis.grid(True, linewidth=0.4, alpha=0.35)
    if equal_aspect:
        axis.set_aspect("equal", adjustable="box")
    figure.savefig(output_path, dpi=200)
    plt.close(figure)


def main():
    args = parse_args()
    path = Path(args.adios_file)
    if not path.exists():
        raise SystemExit(f"ADIOS file does not exist: {path}")

    adios2 = require_module("adios2", "ADIOS2 Python bindings")
    numpy = require_module("numpy", "NumPy")
    matplotlib = require_module("matplotlib", "Matplotlib")
    matplotlib.use("Agg")

    arrays = read_puncture_arrays(adios2, path, numpy)
    puncture_count = len(arrays["y"])
    if puncture_count == 0:
        raise SystemExit(f"No punctures found in {path}")

    yz_path, theta_psi_path = output_paths(args.output_name)
    plot_scatter(
        matplotlib,
        arrays["y"],
        arrays["z"],
        "Y",
        "Z",
        f"Y,Z punctures ({puncture_count} points)",
        yz_path,
        equal_aspect=True,
    )
    plot_scatter(
        matplotlib,
        arrays["theta"],
        arrays["psi"],
        "theta",
        "psi",
        f"theta,psi punctures ({puncture_count} points)",
        theta_psi_path,
    )

    print(f"Wrote {yz_path}")
    print(f"Wrote {theta_psi_path}")


if __name__ == "__main__":
    main()
