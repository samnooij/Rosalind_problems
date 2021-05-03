# Since my algorithm is a bit slow, it may be better to
# move this relatively difficult test case to a separate script.

from pathlib import Path
import sys

sys.path.insert(1, "./")

from scripts.LONG_Genome_assembly_as_shortest_superstring import (
    read_sequences,
    find_overlaps,
)


def test_rosalind_exercise():
    # Make sure that the script produces the right answer as judged by Rosalind
    test_file = "data/rosalind_long5.txt"
    answer_file = "results/long5.txt"

    with open(answer_file, "r") as read_file:
        answer = read_file.readline().strip()

    sequence_dictionary = read_sequences(test_file)
    assert find_overlaps(sequence_dictionary) == answer
