#!/usr/bin/env python3

import csv
from argparse import ArgumentParser
from math import floor
from pathlib import Path
from typing import Iterator, NamedTuple, Optional

import extended_io as eio
import pyhmmer
from xopen import xopen


class HmmHit(NamedTuple):
    score: float
    best_domain_score: float
    included: bool
    sequence: bytes


def msa_proteins_to_digital_seqs(
    seqs: tuple[bytes, ...], alphabet: pyhmmer.easel.Alphabet = pyhmmer.easel.Alphabet.amino()
) -> list[pyhmmer.easel.DigitalSequence]:
    """
    Convert proteins to a list of pyhmmer DigitalSequences.
    :param seqs: byte strings of protein sequences
    :param alphabet: alphabet of the sequences
    :return: list of pyhmmer DigitalSequences
    """
    sequences: list[pyhmmer.easel.DigitalSequence] = []
    for aa in seqs:
        sequences.append(pyhmmer.easel.TextSequence(sequence=aa.decode("UTF-8"), name=aa).digitize(alphabet))
    return sequences


def fasta_proteins_to_digital_seqs(fasta_file: Path) -> list[pyhmmer.easel.DigitalSequence]:
    """
    Convert protein sequences of a FASTA file to a list of digital pyhmmer Sequences.
    :param fasta_file: path to FASTA file
    :return: list of pyhmmer DigitalSequences from the sORF proteins
    """
    sequences: list[pyhmmer.easel.DigitalSequence] = []
    for _, aa in eio.parse_fasta(fasta_file):
        sequences.append(
            pyhmmer.easel.TextSequence(sequence=aa, name=aa.encode("UTF-8")).digitize(pyhmmer.easel.Alphabet.amino())
        )
    return sequences


def perform_hmmsearch(
    proteins: list[pyhmmer.easel.DigitalSequence],
    hmm: pyhmmer.plan7.HMM,
    cutoff: Optional[str] = None,
    threads: int = 1,
) -> list[HmmHit]:
    """
    Perform a hmmsearch on the given proteins with the HMM profile(s).
    :param proteins: list of proteins in pyhmmer DigitalSequence format
    :param hmm: HMM profile(s)
    :param cutoff: cutoff threshold
    :param threads: number of threads
    :return: list of best hits (namedtuple)
    """
    results: list[HmmHit] = []

    search: Iterator[pyhmmer.plan7.TopHits | pyhmmer.plan7.TopHits]
    search = (
        pyhmmer.hmmsearch(hmm, proteins, bit_cutoffs=cutoff, cpus=threads)
        if cutoff
        else pyhmmer.hmmsearch(hmm, proteins, cpus=threads)
    )

    for top_hits in search:
        for hit in top_hits:
            if hit.included:
                results.append(
                    HmmHit(
                        score=hit.score,
                        best_domain_score=hit.best_domain.score,
                        included=hit.included,
                        sequence=hit.name,
                    )
                )
    return results


def calculate_gathering_cutoff(hits) -> tuple[float, float]:
    """
    Calculate the gathering cutoffs from the lowest scoring included sequence.
    :param hits: hits of the proteins against their HMM
    :return: gathering cutoffs (query, domain)
    """
    try:
        gathering_score = min(hit.score for hit in hits if hit.included)
        gathering_domscore = min(hit.best_domain_score for hit in hits if hit.included)
    except ValueError:
        gathering_score, gathering_domscore = 0.0, 0.0
    return floor(gathering_score), floor(gathering_domscore)


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description="Import a SRV abc file and exclude single hits.")
    parser.add_argument("--alignments", "-a", type=Path, help="Input path to directory with alignment files *.aln")
    parser.add_argument("--proteins", "-p", type=Path, help="Input complete sORF FASTA file.")
    parser.add_argument("--output", "-o", type=Path, default=Path("./"), help="Output folder path")
    parser.add_argument("--threads", "-t", type=int, default=1, help="Number of threads")
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    alphabet: pyhmmer.easel.Alphabet = pyhmmer.easel.Alphabet.amino()
    builder: pyhmmer.plan7.Builder = pyhmmer.plan7.Builder(alphabet)
    background: pyhmmer.plan7.Background = pyhmmer.plan7.Background(alphabet)
    hmm: pyhmmer.plan7.HMM

    print("Build HMMs...")
    with (
        open(args.output.joinpath("autoclust.hmm"), "wb") as hmm_file,
        xopen(
            args.output.joinpath("autoclust.clusters..tsv.gz"), "w", compresslevel=9, threads=args.threads
        ) as cluster_file,
    ):
        cluster_tsv = csv.writer(cluster_file, delimiter="\t", lineterminator="\n")
        cluster_tsv.writerow(("protein", "cluster"))

        for aln in eio.find_files(args.alignments, r".*\.afa"):
            with pyhmmer.easel.MSAFile(aln, format="afa") as msa_file:
                msa: pyhmmer.easel.TextMSA = msa_file.read()

                cluster_id: bytes = f"{str(aln).split('/')[-1].lstrip('aln.').rstrip('.afa')}".encode("UTF-8")
                msa.name = cluster_id  # cluster ID: id_subgraph.cluster_size.counter

                digital_msa: pyhmmer.easel.DigitalMSA = msa.digitize(alphabet)
                hmm, _, _ = builder.build_msa(digital_msa, background)

                # field names: http://eddylab.org/software/hmmer/Userguide.pdf page 209
                hmm.accession = cluster_id

                proteins: list[pyhmmer.easel.DigitalSequence] = msa_proteins_to_digital_seqs(
                    digital_msa.names, alphabet=alphabet
                )
                hits: list[HmmHit] = perform_hmmsearch(proteins, hmm, threads=args.threads)

                hmm.cutoffs.gathering = calculate_gathering_cutoff(hits)

                hmm.write(hmm_file)  # , binary=True
                cluster_tsv.writerows(((protein.decode("UTF-8"), cluster_id.decode("UTF-8")) for protein in msa.names))


if __name__ == "__main__":
    main()

# EOF
