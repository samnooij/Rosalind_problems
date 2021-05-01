#! /usr/bin/env python

# A script to find the shortest superstring (contig) for
# a number of DNA sequences provided as fasta file.
# The algorithm works as follows:
# 1. for each sequence, generate 5' and 3' ends of any length
# 2. compare 5' ends to 3' ends of other sequences,
#    and 3' ends to 5' ends of other sequences
# 3. from all overlaps, select the longest,
# 4. merge the sequences with the longest overlap into a new sequence (contig)
# 5. repeat steps 2-4 until no more other sequences exist

import argparse  # Reading command-line arguments (input and output files)
from Bio import SeqIO  # Reading sequences from fasta
from pathlib import Path  # Working with file paths


def parse_arguments():
    """
Parse arguments from the commandline,
that is, input and output file.
    """
    parser = argparse.ArgumentParser(
        description="Takes a fasta file with DNA sequences as input "
        "to return the shortest superstring (contiguous sequence, or contig) "
        "and write that to a text file."
    )

    required = parser.add_argument_group("Required arguments")

    required.add_argument(
        "-i",
        "--input",
        type=str,
        dest="input",
        required=True,
        help="Input file (fasta) with DNA sequences.",
    )

    required.add_argument(
        "-o",
        "--output",
        type=str,
        dest="output",
        required=True,
        help="Output file to which results are written.",
    )

    (args, extra_args) = parser.parse_known_args()

    return args


def read_sequences(input_file):
    """
Given a fasta file as input, read it and return the
sequences as dictionary with IDs as keys and
sequences as values.
(Assumes that each sequence has a unique ID!
Overwrites otherwise.)

    Parameters
----------
input_file : file name
  A `str` with the path to the input file.

Returns
-------
sequence_dictionary : dictionary
  A `dict` with names/IDs as keys and DNA sequences as values.
    """
    sequence_dictionary = {}

    for seq_record in SeqIO.parse(input_file, "fasta"):
        sequence_dictionary[seq_record.id] = seq_record.seq

    return sequence_dictionary


def find_overlaps(seq_dict):
    """
Given a dictionary with DNA sequences as values,
return the shortest contig (supersting).
(Assumptions:
 1. There is at least one way to overlap all sequences,
 2. The overlap is at least one nucleotide from the 5' end 
 and 3' end of two sequences.
 3. The sequences are all in the same orientation (5' to 3')
 )

    Parameters
----------
seq_dict : dictionary with sequences
  A `dict` with DNA sequences as values.

Returns
-------
contig : contiguous sequence
  A `str` with shortest overlapping sequence (contig).
    """

    def list_5_ends(sequence):
        """
        List all possible 5' (left) ends for a given sequence,
        from long to short.
        """
        return [sequence[0:position] for position in range(len(sequence) + 1, 0, -1)]

    def list_3_ends(sequence):
        """
        List all possible 3' (right) ends for a given sequence.
        """
        return [
            sequence[position : len(sequence)] for position in range(0, len(sequence))
        ]

    for sequence in seq_dict.values():
        right_ends = list_5_ends(sequence)

    # This is the part from which to continue the work.
    # I suggest to start by adding the following functions:
    # 1. see if the first sequence is identical to any other:
    #   if so, merge
    # 2. see if the first sequence is identical to any 5' or 3' end:
    #   if so, merge
    # 3. continue with matching the 5' end to all 3' ends of the same length,
    #    and then the 3' end to all other 5' ends:
    #    do this in order from long to short and merge when one has been found
    # experiment with these steps until a satisfactory algorithm has been made.
    #
    # When finished, include test cases that also show how the script works:
    # for example, start by matching sequences 'A' and 'A', then 'AC' and 'CG',
    # add in some interesting but simple cases and slowly increase complexity.
    #
    # Also, when this script works, I would like to implement a version that
    # takes into account reverse complements of sequences as well.
    # (Not particularly useful in this exercise, but may come in handy some time...)


def main():
    """
    Main execution of the script.
    """
    arguments = parse_arguments()

    sequence_dictionary = read_sequences(arguments.input)

    breakpoint()

    return 0


if __name__ == "__main__":
    exit(main())
