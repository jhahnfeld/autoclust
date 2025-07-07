#!/usr/bin/env python3

import os
from argparse import ArgumentParser
from pathlib import Path

import polars as pl


def sequential_melt(df: pl.DataFrame, minsrv: float = 0.0) -> pl.DataFrame:
    """
    Melt the wide Dataframe for sequential for each query to reduce memory usage.
    :param df: Wide Dataframe matrix with subject, query1, query2, ... columns. Values have to SRVs.
    :param minsrv: Report only blast hits with a minimum SRV of minsrv.
    :return: Dataframe in the abc format
    """
    # Create initial dataframe with first query
    abc_df: pl.DataFrame = (
        df.lazy()
        .select(pl.col("subject"), pl.col(df.columns[1]))
        .melt(id_vars="subject", variable_name="query", value_name="srv")
        .drop_nulls()
        .filter(pl.col("srv") >= minsrv)
        .select([pl.col("query"), pl.col("subject"), pl.col("srv")])
        .collect()
    )

    # Create single query melt plans and collect all of them
    lazy_frames = []
    for column in df.columns[2:]:
        lazy_frames.append(
            df.lazy()
            .select(pl.col("subject"), pl.col(column))
            .melt(id_vars="subject", variable_name="query", value_name="srv")
            .drop_nulls()
            .filter(pl.col("srv") >= minsrv)
            .select([pl.col("query"), pl.col("subject"), pl.col("srv")])
        )

    for tmp_df in pl.collect_all(lazy_frames):
        abc_df.vstack(tmp_df, in_place=True)
    return abc_df.rechunk()


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description="Import the diamond  the ")
    parser.add_argument("--input", "-i", type=Path, help="Input a BLAST TSV file with outfmt 6")
    parser.add_argument("--output", "-o", type=Path, default=Path("./"), help="Output folder path")
    parser.add_argument(
        "--mode", "-m", default="blast", type=str, choices=("diamond", "blast"), help='Mode for "diamond" or "blast"'
    )
    parser.add_argument(
        "--minsrv",
        "-n",
        type=int,
        default=0,
        help="Report only blast hits with a minimum SRV of value 30 (= 3.0 * 10).",
    )
    parser.add_argument(
        "--id",
        "-d",
        type=int,
        default=0,
        help="Report only alignments above the given percentage of sequence identity.",
    )
    parser.add_argument(
        "--query_cover",
        "-q",
        type=int,
        default=0,
        help="Report only alignments above the given percentage of query cover.",
    )
    parser.add_argument(
        "--subject_cover",
        "-s",
        type=int,
        default=0,
        help="Report only alignments above the given percentage of subject cover.",
    )
    parser.add_argument("--threads", "-t", type=int, default=1, help="Number of threads")
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    os.makedirs(args.output, exist_ok=True)

    if args.mode not in {"diamond", "blast"}:
        raise Exception('Supported modes: "diamond" or "blast"')

    print(
        f"Settings:\n"
        f"\tmode: {args.mode}\n"
        f"\tid: {args.id}\n"
        f"\tquery_cover: {args.query_cover}\n"
        f"\tsubject_cover: {args.subject_cover}\n"
        f"\tminsrv: {args.minsrv}\n"
    )

    print(f"Load {args.mode} results...")
    if args.mode == "diamond":
        df = (
            pl.read_csv(
                args.input, has_header=False, new_columns=["query", "subject", "bitscore", "evalue"], separator="\t"
            )
            .drop("evalue")
            .group_by("query", "subject")
            .agg(pl.col("bitscore").max())
        )
    else:  # blast
        df = (
            pl.scan_csv(
                args.input,
                has_header=False,
                new_columns=["query", "subject", "identity", "alignment_length", "bitscore", "evalue"],
                schema_overrides={"query": pl.String, "subject": pl.String, "bitscore": pl.Float64},
                separator="\t",
            )
            .filter(
                (pl.col("identity") >= args.id)
                & (pl.col("alignment_length") / pl.col("query").str.len_chars() * 100 >= args.query_cover)
                & (pl.col("alignment_length") / pl.col("subject").str.len_chars() * 100 >= args.subject_cover)
            )
            .select(pl.col("query"), pl.col("subject"), pl.col("bitscore"))
            .group_by("query", "subject")
            .agg(pl.col("bitscore").max())  # Only one hit per sequence (if several hsp hits cause matches)
            .collect()
        )

    print("Pruning...")  # exclude self-hit only
    df = df.filter(
        ~pl.col("query").is_in(
            df.group_by("query")
            .len()
            .filter(pl.col("len") == 1)
            .filter(
                # unique query is in unique subject == self-hit
                pl.col("query").is_in(df.group_by("subject").len().filter(pl.col("len") == 1).get_column("subject"))
            )
            .get_column("query")
        )
    ).rechunk()

    # TODO Replace SRV calculation code, see positive negative filtering code

    print("Transform to matrix...")
    df = df.pivot(index="subject", on="query", values="bitscore")

    print("Calculate and rescale SRVs...")
    df = (
        df.lazy()
        .select(
            pl.col("subject"),
            *(
                pl.col(column) / pl.col(column).max() * 100 for column in df.columns[1:]
            ),  # max per column == query self-hit
        )
        .collect()
    )

    print("Reformat to abc format...")
    df = sequential_melt(df, minsrv=float(args.minsrv))

    print("Pruning & rounding...")  # exclude all self hits
    df = df.filter(pl.col("query") != pl.col("subject")).select(
        pl.col("query"), pl.col("subject"), pl.col("srv").round(4)
    )

    print("Write SRV.abc...")
    df.write_csv(file=str(args.output.joinpath("srv.tsv")), separator="\t", include_header=False)


if __name__ == "__main__":
    main()

# EOF
