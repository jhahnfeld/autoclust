#!/usr/bin/env python3

from argparse import ArgumentParser
from itertools import product
from pathlib import Path

import polars as pl


def read_mcl_cls(file: Path) -> dict[str, int]:
    cluster_map: dict[str, int] = dict()
    with open(file, "r") as infile:
        for i, line in enumerate(infile):
            proteins: list[str] = line.strip().split()
            for protein, cls in product(proteins, [str(i)]):
                cluster_map[protein] = int(cls)
    return cluster_map


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(
        description="Import a mcl abc file and tab files with clusters and export them for use in cytoscape."
    )
    parser.add_argument("--abc", "-a", type=Path, help="Input mcl tab file with clusters")
    parser.add_argument(
        "--autoclust", "-u", type=Path, help="MCL clusters for inflation values chosen for binned subgraphs."
    )
    # parser.add_argument(
    #     "--cluster",
    #     "-c",
    #     nargs="+",
    #     type=list[list[str]],
    #     dest="cluster",
    #     help="MCL clusters for different inflation values.",
    # )
    parser.add_argument("--output", "-o", type=Path, default=Path("./"), help="Output folder path")
    parser.add_argument("--threads", "-t", type=int, default=1, help="Number of threads")
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    graph_df = pl.read_csv(args.abc, has_header=False, new_columns=["query", "subject", "srv"], separator="\t")

    # inflations: list[str] = ["q-autoclust"]
    clusters = read_mcl_cls(args.autoclust)
    graph_df = graph_df.with_columns(pl.col("query").replace(clusters).alias("q-autoclust"))

    # for clustering in args.cluster:
    #     clustering: Path = Path("".join(clustering))
    #     inflation: str = ".".join(str(clustering).split("/")[-1].split(".")[-3:-1])
    #     inflations.append(f"q-{inflation}")
    #     clusters = read_mcl_cls(clustering)
    #
    #     graph_df = graph_df.with_columns(pl.col("query").replace(clusters).alias(f"q-{inflation}"))

    graph_df.write_csv(args.output.joinpath("srv.clustered.tsv"), separator="\t")


if __name__ == "__main__":
    main()

# EOF
