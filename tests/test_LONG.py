from pathlib import Path
import sys

sys.path.insert(1, "./")

from scripts.LONG_Genome_assembly_as_shortest_superstring import (
    read_sequences,
    find_overlaps,
)

sample_dictionary = {
    "Rosalind_56": "ATTAGACCTG",
    "Rosalind_57": "CCTGCCGGAA",
    "Rosalind_58": "AGACCTGCCG",
    "Rosalind_59": "GCCGGAATAC",
}


def test_read_sequences():
    # Make sure the sample dataset works as intended
    sequence_dictionary = read_sequences(
        "data/Example_genome_assembly_as_shortest_superstring.fasta"
    )

    assert len(sequence_dictionary) == 4
    assert sequence_dictionary == sample_dictionary


def test_find_overlaps():
    # Make sure the sample dataset works as intended
    assert find_overlaps(sample_dictionary) == "ATTAGACCTGCCGGAATAC"


# Now start a series of simple tests


def test_single_sequence():
    # This should return the sequence itself
    assert find_overlaps({"single_sequence": "ACGT"}) == "ACGT"


def test_simple_deduplication():
    # This should return on of the two identical sequences
    simple_duplicates = {"one": "A", "two": "A"}
    assert find_overlaps(simple_duplicates) == "A"


def test_two_simple_sequences():
    two_simple_sequences = {"one": "AC", "two": "CG"}
    assert find_overlaps(two_simple_sequences) == "ACG"


def test_simple_sequences_with_duplicates():
    # This should return a single instance of the longest supersequence
    simple_with_duplicates = {
        "one": "AC",
        "two": "CG",
        "three": "GT",
        "one-dup": "AC",
        "two-dup": "CG",
        "three-dup": "GT",
    }
    assert find_overlaps(simple_with_duplicates) == "ACGT"


def test_longest_overlap():
    # This should return the longest possible overlap =
    # shortest possible subsequence between sequences.
    longest_overlap1 = {"one": "ACCCC", "two": "CCCCG"}
    assert find_overlaps(longest_overlap1) == "ACCCCG"

    longest_overlap2 = {"three": "GGGGT", "four": "AGGGG"}
    assert find_overlaps(longest_overlap2) == "AGGGGT"
    # Matching should work in both directions


# Since there is not yet support for sequences that do
# not overlap, there are no tests with only different
# sequences, for instance "A" and "C".
# This would currently result in an infinite loop.
# def test_non_overlapping_sequences():
#    no_overlap = {"one": "A", "two": "C"}
#    assert find_overlaps(no_overlap) == no_overlap #??
