#!/usr/bin/env python3

import argparse
import shlex


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate calc_dxdz_dy.py commands for a timestep range."
    )
    parser.add_argument("apar_file")
    parser.add_argument("grid_file")
    parser.add_argument("n0", type=int)
    parser.add_argument("n1", type=int)
    parser.add_argument("outputname")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.n0 < 0:
        raise SystemExit("n0 must be non-negative")
    if args.n1 < args.n0:
        raise SystemExit("n1 must be greater than or equal to n0")

    for timestep in range(args.n0, args.n1 + 1):
        output_file = f"{args.outputname}.{timestep:05d}.nc"
        command = [
            "python3",
            "calc_dxdz_dy.py",
            args.grid_file,
            args.apar_file,
            "1",
            str(timestep),
            output_file,
        ]
        print(" ".join(shlex.quote(item) for item in command))


if __name__ == "__main__":
    main()
