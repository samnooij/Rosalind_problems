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


class font:
    """
Write print messages with bold and/or coloured letters.
Reset to normal with 'font.reset'.
Start print statement with e.g. 'font.bold, font.red'
to make bold red letters.
    """

    reset = "\033[0m"
    bold = "\033[01m"
    red = "\033[31m"
    # See more examples at
    # https://www.geeksforgeeks.org/print-colors-python-terminal/


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

    def highlight_substring(shorter, longer):
        """
        Print a shorter sequence that is a subsequence
        of a longer sequence, and highlight the shorter
        sequence within the longer sequence in bold red.
        """
        length = len(shorter)
        start = longer.index(shorter)
        end = start + length

        print("%s\nis a subsequence of" % shorter)
        print(
            longer[:start]
            + font.bold
            + font.red
            + longer[start:end]
            + font.reset
            + longer[end:]
            + "\n"
        )

    def highlight_overlap(left, right, length):
        """
        Highlight the overlapping part between a left and
        right sequence in bold red print.
        """
        left_start = len(left) - length

        print(
            "Sequences:\n"
            + left[:left_start]
            + font.bold
            + font.red
            + left[left_start:]
            + font.reset
            + " and "
            + font.bold
            + font.red
            + right[:length]
            + font.reset
            + right[length:]
            + " overlap!"
        )

    def list_left_ends(sequence):
        """
        List all possible left (probably 5') ends for a given sequence,
        from long to short.
        """
        return [sequence[0:position] for position in range(len(sequence) + 1, 0, -1)]

    def list_right_ends(sequence):
        """
        List all possible right (probably 3') ends for a given sequence.
        """
        return [
            sequence[position : len(sequence)] for position in range(0, len(sequence))
        ]

    seq_dict_copy = seq_dict.copy()

    current_id = list(seq_dict_copy.keys())[0]
    current_sequence = seq_dict_copy.pop(current_id)

    # First, check if there are any other sequences,
    if len(seq_dict_copy) == 0:
        # If there are none, return the current sequence
        return current_sequence

    # then check if there is a duplicate
    elif current_sequence in seq_dict_copy.values():
        # And if there is, restart the process
        # excluding the duplicate 'current_sequence'
        print("There is a duplicate of %s" % (current_sequence))

        return find_overlaps(seq_dict_copy)

    else:
        # If there is no duplicate, continue by looking
        # for duplicates as substrings: is the current sequence
        # a substring of a longer sequence in the dictionary,
        # or is any of the sequences in the dictionary a
        # subsequence of the current sequence?
        for name, sequence in seq_dict.items():
            if len(sequence) > len(current_sequence):
                # If the sequence from the dictionary is longer,
                # current_sequence may be a substring
                if current_sequence in sequence:
                    highlight_substring(current_sequence, sequence)

                    return find_overlaps(seq_dict_copy)
                    # Continue without current sequence

            elif len(sequence) < len(current_sequence):
                # Or if the current sequence is longer,
                # the sequence from the dictionary may be
                # a substring of the current.
                if sequence in current_sequence:
                    highlight_substring(sequence, current_sequence)

                    del seq_dict_copy[name]  # Remove the shorter sequence
                    seq_dict_copy.update({name + "-deduplicated": current_sequence})
                    return find_overlaps(seq_dict_copy)

            else:
                # If neither is longer, neither can be a substring
                # of the other.
                pass

        # Finally, look for overlaps in both left and right (
        # 5' and 3') ends.

        current_left_ends = list_left_ends(current_sequence)
        current_right_ends = list_right_ends(current_sequence)

        for index in range(1, len(current_sequence)):
            left_end = current_left_ends[index]
            right_end = current_right_ends[index]
            # For each left end and right end,

            for name, sequence in seq_dict_copy.items():
                if sequence.endswith(left_end):
                    # If another sequence ends with the current
                    # sequence's left end, there is an overlap
                    highlight_overlap(sequence, current_sequence, len(left_end))
                    merged_sequence = current_sequence[len(left_end) :] + sequence

                    print("New sequence: %s\n" % merged_sequence)
                    seq_dict_copy.update({name + "-merged": merged_sequence})
                    return find_overlaps(seq_dict_copy)
                    # It is assumed there is only one longest match
                    # in the dataset, so stop looking further.

                elif sequence.startswith(right_end):
                    # And if another sequence ends with the
                    # current sequence's right end, there is an overap
                    highlight_overlap(current_sequence, sequence, len(right_end))
                    merged_sequence = current_sequence[: -len(right_end)] + sequence

                    print("New sequence: %s\n" % merged_sequence)
                    seq_dict_copy.update({name + "-merged": merged_sequence})
                    return find_overlaps(seq_dict_copy)
                    # It is assumed there is only one longest match
                    # in the dataset, so stop looking further.

                else:
                    # If the sequence is not a duplicate and has no overlap,
                    # throw it out (for now):
                    # print(sorted(list(seq_dict_copy.keys())))
                    # del seq_dict_copy[name]
                    merged_sequence = current_sequence

    seq_dict_copy.update({"merged_sequence": merged_sequence})

    return find_overlaps(seq_dict_copy)

    # When finished, include test cases that also show how the script works:
    # for example, start by matching sequences 'A' and 'A', then 'AC' and 'CG',
    # add in some interesting but simple cases and slowly increase complexity.
    #
    # Also, when this script works, I would like to implement a version that
    # takes into account reverse complements of sequences as well.
    # (Not particularly useful in this exercise, but may come in handy some time...)
    # After that, it may be nice to add a version that does not require
    # all sequences to have at least some overlap and just returns multiple
    # sequences if some sequences share no overlapping parts.
    # (This should probably also use a cut-off for overlap length:
    # overlaps of length 1 are not that hard to find, but are often meaningless
    # when applied to actual genomics data.)


def main():
    """
    Main execution of the script.
    """
    arguments = parse_arguments()

    sequence_dictionary = read_sequences(arguments.input)

    shortest_supersequence = find_overlaps(sequence_dictionary)

    print(
        "-----\nThe shortest supersequence is:\n"
        + font.bold
        + shortest_supersequence
        + font.reset
    )

    with open(arguments.output, "w") as results_file:
        # For a fasta file, add an identifier line
        if arguments.output.endswith(".fasta"):
            results_file.write(">Assembled_sequence\n%s\n" % shortest_supersequence)

        else:
            # Otherwise, write only the sequence as plain text
            results_file.write("%s\n" % shortest_supersequence)

    return 0


if __name__ == "__main__":
    exit(main())

## Minimum required features:
# 1. Compare ends of sequences, long to short
#  --> merge to first = longest match, there should be only one
# 2.
#
## Bonus features for extended usability:
# 1. Collapse all duplicates
#  --> include duplicates that are subsequences of longer sequences
# 2. Continue searching after an end has been matched: there may be multiple options
#  --> in case of multiple options, try all or stop extending?
# 3. Include reverse complements as option
