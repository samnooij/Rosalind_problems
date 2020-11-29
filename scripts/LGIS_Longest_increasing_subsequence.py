#! /usr/bin/env python

import argparse
from bisect import bisect_left, bisect_right
from functools import cmp_to_key


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


# Function 'longest_subsequence' kindly provided by arekolek on stackoverflow:
# https://stackoverflow.com/a/38337443
def longest_subsequence(
    seq, mode="strictly", order="increasing", key=None, index=False
):
    """
Return the longest increasing subsequence of `seq`.

Parameters
----------
seq : sequence object
  Can be any sequence, like `str`, `list`, `numpy.array`.
mode : {'strict', 'strictly', 'weak', 'weakly'}, optional
  If set to 'strict', the subsequence will contain unique elements.
  Using 'weak' an element can be repeated many times.
  Modes ending in -ly serve as a convenience to use with `order` parameter,
  because `longest_sequence(seq, 'weakly', 'increasing')` reads better.
  The default is 'strict'.
order : {'increasing', 'decreasing'}, optional
  By default return the longest increasing subsequence, but it is possible
  to return the longest decreasing sequence as well.
key : function, optional
  Specifies a function of one argument that is used to extract a comparison
  key from each list element (e.g., `str.lower`, `lambda x: x[0]`).
  The default value is `None` (compare the elements directly).
index : bool, optional
  If set to `True`, return the indices of the subsequence, otherwise return
  the elements. Default is `False`.

Returns
-------
elements : list, optional
  A `list` of elements of the longest subsequence.
  Returned by default and when `index` is set to `False`.
indices : list, optional
  A `list` of indices pointing to elements in the longest subsequence.
  Returned when `index` is set to `True`.
    """
    bisect = bisect_left if mode.startswith("strict") else bisect_right

    # compute keys for comparison just once
    rank = seq if key is None else map(key, seq)
    if order == "decreasing":
        rank = map(cmp_to_key(lambda x, y: 1 if x < y else 0 if x == y else -1), rank)
    rank = list(rank)

    if not rank:
        return []

    lastoflength = [0]  # end position of subsequence with given length
    predecessor = [None]  # penultimate element of l.i.s. ending at given position

    for i in range(1, len(seq)):
        # seq[i] can extend a subsequence that ends with a lesser (or equal) element
        j = bisect([rank[k] for k in lastoflength], rank[i])
        # update existing subsequence of length j or extend the longest
        try:
            lastoflength[j] = i
        except:
            lastoflength.append(i)
        # remember element before seq[i] in the subsequence
        predecessor.append(lastoflength[j - 1] if j > 0 else None)

    # trace indices [p^n(i), ..., p(p(i)), p(i), i], where n=len(lastoflength)-1
    def trace(i):
        if i is not None:
            yield from trace(predecessor[i])
            yield i

    indices = trace(lastoflength[-1])

    return list(indices) if index else [seq[i] for i in indices]


def main():
    """
    Main execution of the script.
    """
    arguments = parse_arguments()
    (length, sequence) = read_length_and_sequence(arguments.input)

    increasing_subsequence = longest_subsequence(
        seq=sequence, mode="strictly", order="increasing"
    )
    decreasing_subsequence = longest_subsequence(
        seq=sequence, mode="strictly", order="decreasing"
    )

    return 0


if __name__ == "__main__":
    exit(main())
