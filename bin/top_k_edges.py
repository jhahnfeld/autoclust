#!/usr/bin/env python3

import os
import subprocess as sp
import sys
from argparse import ArgumentParser
from pathlib import Path


def mcx_query(imx: Path, out_path: Path, top_k: int = 1000, low_k: int = 0, step: int = 100, threads: int = 1):
    cmd = [
        "mcx",
        "query",
        "-imx",
        str(imx),
        "-vary-knn",
        f"{low_k}/{top_k}/{step}",
        "-t",
        str(threads)
    ]
    with open(out_path, mode='w') as f:
        proc = sp.run(cmd, env=os.environ.copy(), stdout=f, universal_newlines=True)  #, stdout=sp.PIPE, stderr=sp.PIPE
        if proc.returncode != 0:
            sys.stderr.write(f"Command:\n{' '.join(cmd)}")
            sys.stderr.write(f"stdout={proc.stdout}, stderr={proc.stderr}")
            raise Exception(f"mcx error! error code: {proc.returncode}")


def get_new_k(mcx_path: Path) -> int:
    start: int = 0
    last_k: int = 1000
    with open(mcx_path, mode="r") as fh:
        for line in fh:
            splitted_line = line.strip().split()
            if (
                splitted_line[0]
                == "----------------------------------------------------------------------------------------------"
            ):
                start += 1
            elif (
                splitted_line[0]
                == "---------------------------------------------------------------------------------------------"
            ):
                continue
            elif start == 2:
                current_k: int = int(splitted_line[14])
                num_singletons: int = int(splitted_line[3])
                if num_singletons >= 1:
                    return last_k
                else:
                    last_k = current_k
        else:
            return last_k


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(
        description="Import a mcl abc file and tab files with clusters and export them for use in cytoscape."
    )
    parser.add_argument("--imx", "-i", type=Path, default=Path("./"), help="Input mcl tab file with clusters")
    parser.add_argument("--output", "-o", type=Path, default=Path("./"), help="Output folder path")
    parser.add_argument("--threads", "-t", type=int, default=1, help="Number of threads")
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    out_path: Path = args.output.joinpath("mille.txt")
    mcx_query(args.imx, out_path, 1000, 0, 100, args.threads)
    mille_k: int = get_new_k(out_path)

    out_path: Path = args.output.joinpath("centum.txt")
    mcx_query(args.imx, out_path, mille_k, mille_k - 100, 10, args.threads)
    centum_k: int = get_new_k(out_path)

    out_path: Path = args.output.joinpath("decem.txt")
    mcx_query(args.imx, out_path, centum_k, centum_k - 10, 1, args.threads)
    decem_k: int = get_new_k(out_path)

    # TODO if decem ist too low print a higher decem instead loook up mcl doku
    # If your input graph is extremely dense, with an average node degree (i.e. the number of neighbours per node) that
    # is somewhere above 500, you may need to filter the input graph by removing edges, for example by using one of -tf '#ceilnb()' or -tf '#knn()'.
    print(decem_k)


if __name__ == "__main__":
    main()
