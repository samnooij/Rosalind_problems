from pathlib import Path
import sys

sys.path.insert(1, "./")

from scripts.LGIS_Longest_increasing_subsequence import (
    read_length_and_sequence,
    longest_subsequence,
    write_lists_as_separate_lines,
)

# function import method thanks to addicted on stackoverflow:
# https://stackoverflow.com/questions/4383571/importing-files-from-different-folder#comment104909466_40612922

test_sequence = ["5", "1", "4", "2", "3"]


def test_read_length_and_sequence():
    with open("data/Example_longest_increasing_subsequence.txt", "r") as read_file:
        length = int(read_file.readline())
        sequence = list(read_file.readline().split())

    assert length == 5
    assert sequence == test_sequence


def test_longest_subsequence():
    assert longest_subsequence(
        seq=test_sequence, mode="strictly", order="increasing"
    ) == ["1", "2", "3"]
    assert longest_subsequence(
        seq=test_sequence, mode="strictly", order="decreasing"
    ) == ["5", "4", "3"] or longest_subsequence(
        seq=test_sequence, mode="strictly", order="decreasing"
    ) == [
        "5",
        "4",
        "2",
    ]
    # In this example, there are two possibilities for the longest decreasing subsequence.


def test_write_lists_as_separate_lines():
    ## First test:
    # This test takes into account lists that consist of integers
    # and characters. It tests two lists.

    write_lists_as_separate_lines(
        lists=[[1, 2, 3], ["Testing", "my", "function"],],
        output_file="test_write_lists1.txt",
    )

    # Make sure the output file has been created
    assert Path("test_write_lists1.txt").exists()

    # Check the content of the file
    with open("test_write_lists1.txt", "r") as read_file:
        first_list = read_file.readline().strip()
        second_list = read_file.readline().strip()

        assert first_list == "1 2 3"
        assert second_list == "Testing my function"

    # Remove the file after checking
    Path("test_write_lists1.txt").unlink()
