#! /usr/bin/env python

import argparse


def parse_arguments():
    """
Parse arguments from the commandline,
that is, input and output file.
    """
    parser = argparse.ArgumentParser(description="Read input parameters.")

    required = parser.add_argument_group("Required arguments")

    required.add_argument(
        "-i",
        "--input",
        type=str,
        dest="input",
        required=True,
        help="Input file with lenth and sequence of integers.",
    )

    (args, extra_args) = parser.parse_known_args()

    return args


def read_length_and_sequence(input_file):
    """
Given an input file with a length (integer) on the first line,
and a sequence of integers on the second line,
read the length and sequence as variables (integer and list).
    """
    with open(input_file, "r") as read_file:
        length = int(read_file.readline())
        sequence = list(read_file.readline().split())

    return (length, sequence)


def main():
    """
    Main execution of the script.
    """
    arguments = parse_arguments()
    (length, sequence) = read_length_and_sequence(arguments.input)

    return 0


if __name__ == "__main__":
    exit(main())
