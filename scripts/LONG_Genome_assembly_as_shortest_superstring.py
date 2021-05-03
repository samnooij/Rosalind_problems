#! /usr/bin/env python3

# Requires Python 3.7+ to work with ordered dictionaries!
# https://docs.python.org/3/library/stdtypes.html#mapping-types-dict

# A script to find the shortest superstring (contig) for
# a number of DNA sequences provided as fasta file.
# The algorithm works as follows (pseudocode):
#
#    - read sequences from a fasta file, store as dictionary
#  (1) take out the first element of the dictionary
#  (2) if there are no other sequences (there is only the first one):
#        return the sequence
#  (3) else, if the sequence is identical to any other sequence:
#        discard the sequence, start from (1)
#  (4) else, if any other sequence is a substring of the current sequence,
#      or the current sequence is a substring of any other sequence:
#        discard the shorter one, start from (1)
#      else, there are no complete duplicates
#  (5) for each left or right end of the sequence (from long to short),
#      (6) check with each other sequence:
#           (7) if the ends are at least half the sequences' lengths
#                (8a) if the left end matches:
#                        merge sequences, remove the other sequence, start from (1)
#                (8b)  if the right end matches:
#                        merge the sequences, remove the other sequence, start from (1)
#                    if neither end matches:
#                        continue with next (other) sequence
#                if the overlap is not at least half as long as the sequences:
#                    continue until no left and right ends are left
#            (9) if no overlaps can be found:
#                    put sequence in end of line, start from (1)
#
# Numbered steps are labelled in the script as comments.

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
return the shortest contig (superstring).
(Assumptions:
 1. There is at least one way to overlap all sequences,
 2. The overlap is at half of either sequence's length.
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
        List all possible right (probably 3') ends for a given sequence,
        from long to short.
        """
        return [
            sequence[position : len(sequence)] for position in range(0, len(sequence))
        ]

    seq_dict_copy = seq_dict.copy()

    # (1) Take the first element:
    current_id = list(seq_dict_copy)[0]
    current_sequence = seq_dict_copy.pop(current_id)

    # (2) If there are no other sequences:
    if len(seq_dict_copy) == 0:
        # Return the current sequence
        return current_sequence

    # (3) If the sequence is identical to any other:
    elif current_sequence in seq_dict_copy.values():
        print("There is a duplicate of %s" % (current_sequence))

        # Start from (1), discarding the current sequence
        return find_overlaps(seq_dict_copy)

    else:
        # (4)
        # If there is no duplicate, continue by looking
        # for duplicates as substrings: is the current sequence
        # a substring of a longer sequence in the dictionary,
        # or is any of the sequences in the dictionary a
        # subsequence of the current sequence?
        for name, sequence in seq_dict_copy.items():
            if len(sequence) > len(current_sequence):
                # If the sequence from the dictionary is longer,
                # current_sequence may be a substring
                if current_sequence in sequence:
                    highlight_substring(current_sequence, sequence)

                    return find_overlaps(seq_dict_copy)
                    # Start from (1) without current sequence

            elif len(sequence) < len(current_sequence):
                # Or if the current sequence is longer,
                # the sequence from the dictionary may be
                # a substring of the current.
                if sequence in current_sequence:
                    highlight_substring(sequence, current_sequence)

                    del seq_dict_copy[name]  # Remove the shorter sequence

                    new_dict = {name + "-deduplicated": current_sequence}
                    new_dict.update(seq_dict_copy)
                    # Add the current sequence at the start of a new dictionary

                    return find_overlaps(new_dict)
                    # Start from (1) without the shorter duplicate

            else:
                # If neither is longer, neither can be a substring
                # of the other.
                pass

        # Finally, look for overlaps in both left and right (
        # 5' and 3') ends.

        current_left_ends = list_left_ends(current_sequence)
        current_right_ends = list_right_ends(current_sequence)

        # (5) For each left and right end of the sequence
        # (from long to short).
        for index in range(1, len(current_sequence)):
            left_end = current_left_ends[index]
            right_end = current_right_ends[index]

            # (6) Compare to each other sequence:
            for name, sequence in seq_dict_copy.items():
                # (7) If the overlap part is at least
                # half as long as either of the sequences:
                if (
                    len(left_end) >= len(sequence) / 2
                    or len(left_end) >= len(current_sequence) / 2
                ):
                    # Left and right should be equally long, so check only one.
                    # If the overlap is at least half as long as both sequences,
                    # continue checking:

                    # (8a) If the left end matches the other sequence:
                    if sequence.endswith(left_end):
                        highlight_overlap(sequence, current_sequence, len(left_end))
                        merged_sequence = sequence + current_sequence[len(left_end) :]

                        print("New sequence: %s\n" % merged_sequence)

                        new_dict = {name + "-merged": merged_sequence}
                        new_dict.update(seq_dict_copy)
                        # Add the merged sequence at the start of a new dictionary

                        del new_dict[name]  # Remove the matched sequence

                        return find_overlaps(new_dict)
                        # Start from (1) with the two sequences merged.

                    # (8b) If the right end matches the other sequence:
                    elif sequence.startswith(right_end):
                        highlight_overlap(current_sequence, sequence, len(right_end))
                        merged_sequence = current_sequence[: -len(right_end)] + sequence

                        print("New sequence: %s\n" % merged_sequence)

                        new_dict = {name + "-merged": merged_sequence}
                        new_dict.update(seq_dict_copy)
                        # Add the merged sequence at the start of a new dictionary

                        del new_dict[name]  # Remove the matched sequence

                        return find_overlaps(new_dict)
                        # Start from (1) with the two sequences merged.

                    else:
                        # If the sequence is not a duplicate and has no overlap,
                        # continue with the next sequence.
                        pass

                else:
                    # If the overlap part is shorter than half a sequence:
                    pass  # Do not consider these right/left ends further.

        # (9) If no overlaps are found:
        print(
            "There appears to be no overlap between %s and any other sequence."
            % current_sequence
        )

        # Put the sequence at the end of the dictionary...
        print("Retrying with next sequence...")
        seq_dict_copy.update({current_id: current_sequence})

    # ...and start from (1)
    return find_overlaps(seq_dict_copy)


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
