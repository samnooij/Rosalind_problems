#!/usr/bin/env python
# coding: utf-8

# # Inferring mRNA from protein, from [Rosalind.info](https://www.rosalind.info)
# 
# (Text copied from http://rosalind.info/problems/mrna/)
# 
# <div class="problem-statement problem-statement-bordered" problem="130">
#     <blockquote>
# <h2 id="pitfalls-of-reversing-translation">Pitfalls of Reversing Translation</h2>
# <p>When researchers discover a new <a class="term" href="/glossary/protein/" id="term-206" rel="tooltip" title="The functional unit of the cell.">protein</a>, they would like to infer
# the strand of <a class="term" href="/glossary/messenger-rna/" id="term-266" rel="tooltip" title="
# An RNA molecule that serves as the blueprint for translation into protein.">mRNA</a> from which this protein could have been <a class="term" href="/glossary/translation/" id="term-247" rel="tooltip" title="
# The process by which mRNA is converted into a peptide chain for the creation of a protein.">translated</a>,
# thus allowing them to locate <a class="term" href="/glossary/gene/" id="term-363" rel="tooltip" title="
# An interval of DNA whose nucleotides are translated into a polypeptide for protein creation.">genes</a> associated with this protein on the <a class="term" href="/glossary/genome/" id="term-360" rel="tooltip" title="
# The collection of all of an organism's DNA taken from all its chromosomes.">genome</a>.</p>
# <p>Unfortunately, although any <a class="term" href="/glossary/rna-string/" id="term-350" rel="tooltip" title="
# A string constructed from the alphabet {A, C, G, U}.">RNA string</a> can be translated into a unique <a class="term" href="/glossary/protein-string/" id="term-371" rel="tooltip" title="
# A string composed of symbols taken from the English alphabet less B, J, O, U, X, and Z;
# representing a peptide chain formed from amino acids.">protein string</a>,
# reversing the process yields a huge number of possible RNA strings from a single
# protein string because most amino acids correspond to multiple RNA <a class="term" href="/glossary/codon/" id="term-259" rel="tooltip" title="
# A triplet of contiguous nucleotides.">codons</a>
# (see the <a class="term" href="/glossary/rna-codon-table/" id="term-217" rel="tooltip" title="A table indicating the translation of individual RNA codons into amino acids for the
# purpose of protein creation.">RNA Codon Table</a>).</p>
# <p>Because of memory considerations, most data formats that are built into languages have
# upper bounds on how large an integer can be: in some versions of Python,
# an "int" variable may be required to be no larger than <mathjax>$2^{31} -1$</mathjax>, or 2,147,483,647.
# As a result, to deal with very large numbers in Rosalind, we need to devise a system
# that allows us to manipulate large numbers without actually having to store large numbers.</p>
# </blockquote>
# <h2 id="problem">Problem</h2>
# <p>For positive integers <mathjax>$a$</mathjax> and <mathjax>$n$</mathjax>, <mathjax>$a$</mathjax> <a class="term new" href="/glossary/modular-arithmetic/" id="term-577" rel="tooltip" title="New term: 
# The study of arithmetic on integer remainders.">modulo</a> <mathjax>$n$</mathjax> (written <mathjax>$a\mod n$</mathjax> in shorthand) is the remainder
# when <mathjax>$a$</mathjax> is divided by <mathjax>$n$</mathjax>.  For example, <mathjax>$29 \mod 11 = 7$</mathjax> because <mathjax>$29 = 11 \times 2 + 7$</mathjax>.</p>
# <p><a class="term new" href="/glossary/modular-arithmetic/" id="term-577" rel="tooltip" title="New term: 
# The study of arithmetic on integer remainders.">Modular arithmetic</a> is the study of addition, subtraction, multiplication, and division
# with respect to the modulo operation.  We say that <mathjax>$a$</mathjax> and <mathjax>$b$</mathjax> are <a class="term new" href="/glossary/modular-arithmetic/" id="term-577" rel="tooltip" title="New term: 
# The study of arithmetic on integer remainders.">congruent</a> modulo <mathjax>$n$</mathjax>
# if <mathjax>$a \mod n = b \mod n$</mathjax>; in this case, we use the notation <mathjax>$a \equiv b \mod n$</mathjax>.</p>
# <p>Two useful facts in modular arithmetic are that if <mathjax>$a \equiv b \mod n$</mathjax> and <mathjax>$c \equiv d \mod n$</mathjax>,
# then <mathjax>$a+c \equiv b+d \mod n$</mathjax> and <mathjax>$a \times c \equiv b \times d \mod n$</mathjax>.  To check your understanding of these rules,
# you may wish to verify these relationships for <mathjax>$a = 29$</mathjax>, <mathjax>$b = 73$</mathjax>, <mathjax>$c = 10$</mathjax>, <mathjax>$d = 32$</mathjax>, and <mathjax>$n = 11$</mathjax>.</p>
# <p>As you will see in this exercise, some Rosalind problems will ask for a (very large)
# integer solution modulo a smaller number to avoid the computational pitfalls that arise with
# storing such large numbers.</p>
# <p><span class="given-return">Given:</span> A <a class="term" href="/glossary/protein-string/" id="term-371" rel="tooltip" title="
# A string composed of symbols taken from the English alphabet less B, J, O, U, X, and Z;
# representing a peptide chain formed from amino acids.">protein string</a> of length at most 1000 <a class="term" href="/glossary/amino-acid/" id="term-198" rel="tooltip" title="
# The monomer unit for proteins; the same 20 amino acids commonly occur in most species.">aa</a>.</p>
# <p><span class="given-return">Return:</span> The total number of different RNA strings from which the protein could have been translated, modulo 1,000,000.
# (Don't neglect the importance of the <a class="term" href="/glossary/stop-codon/" id="term-261" rel="tooltip" title="
# One of three possible RNA codons that indicate the termination of protein translation.">stop codon</a> in protein translation.)</p>
# <h2 id="sample-dataset">Sample Dataset</h2>
# <div class="codehilite"><pre>MA
# </pre></div>
# 
# 
# <h2 id="sample-output">Sample Output</h2>
# <div class="codehilite"><pre>12
# </pre></div>
# 
# 
# <blockquote>
# <h2 id="hint">Hint</h2>
# <p>What does it mean intuitively to take a number modulo 1,000,000?</p>
# </blockquote>
# </div>

# ## My interpretation/reasoning
# 
# 1. First the hint: Intuitively, to me this means that the final number is smaller than the modulo number. So a (big) number modulo 1,000,000 will return a number smaller than 1,000,000, which is safely below the integer limit of $2^{31}-1$.
# 
# 2. This problem is about 'backtranslating' amino acid sequences (proteins) into nucleotide sequences (mRNA).
# 
# 3. However, the exercise does not ask to make actual mRNA sequences.
# 
# 4. Instead, I get to calculate how many options exist for the backtranslation.
# 
# 5. For some amino acid there is only one codon, for some there are multiple. (Note that there are also multiple stop codons! This is not written in the amino acid sequence, but does add 3 possibilities.)
# 
# 6. The final answer should be the product of the possibilities per amino acid, e.g. "MA" reads: 1 codon for M * 4 codons for A * 3 stop codons = 1 * 4 * 3 = 12 possible mRNA sequences.
# 
# _Note: Rosalind has included this RNA codon table as reminder:_
# 
# ```
# UUU F      CUU L      AUU I      GUU V
# UUC F      CUC L      AUC I      GUC V
# UUA L      CUA L      AUA I      GUA V
# UUG L      CUG L      AUG M      GUG V
# UCU S      CCU P      ACU T      GCU A
# UCC S      CCC P      ACC T      GCC A
# UCA S      CCA P      ACA T      GCA A
# UCG S      CCG P      ACG T      GCG A
# UAU Y      CAU H      AAU N      GAU D
# UAC Y      CAC H      AAC N      GAC D
# UAA Stop   CAA Q      AAA K      GAA E
# UAG Stop   CAG Q      AAG K      GAG E
# UGU C      CGU R      AGU S      GGU G
# UGC C      CGC R      AGC S      GGC G
# UGA Stop   CGA R      AGA R      GGA G
# UGG W      CGG R      AGG R      GGG G
# ```
# 
# To put this into a script, I want to have:
#  - A list of codons per amino acid
#  - A function to read the amino acid sequence (from a text file)
#  - Then for each amino acid in the sequence, I want to multiply the number of possibilities by the corresponding number of possibilities (starting with 1)
#  - The final product should be _first_ multiplied by 3 to add in the stop codon, _then_ take the modulo 1,000,000 to produce the final answer.
#  
# Sounds not too difficult, right? Let's hope it work out!

# In[11]:


def read_protein_sequence(input_file):
    """
Read a protein sequence from an input file.
Return the sequence as string (text).

Note: I am assuming exactly 1 line of text.
    """
    with open(input_file, 'r') as read_file:
        for line in read_file:
            sequence = line.strip()
            
    return(sequence)


# In[2]:


test_sequence = read_protein_sequence("data/Example_Inferring_mRNA_from_protein.txt")

print(test_sequence)


# In[4]:


number_of_codons = {
    "A": 4,
    "C": 2,
    "D": 2,
    "E": 2,
    "F": 2,
    "G": 4,
    "H": 2,
    "I": 3,
    "K": 2,
    "L": 6,
    "M": 1,
    "N": 2,
    "P": 4,
    "Q": 2,
    "R": 6,
    "S": 6,
    "T": 4,
    "V": 4,
    "W": 1,
    "Y": 2,
    "Stop": 3
}

#As a little check: the sum of these numbers should be 4 * 16 = 64:

print(sum(number_of_codons.values()))


# In[5]:


from functools import reduce
from operator import mul
#Use these functions to multiple the contents of a list


# In[6]:


def calculate_mrna_possibilities(aa_sequence):
    """
Given a protein sequence (amino acids),
calculate the number of mRNA (nucleotide) sequences
that can translate into that protein sequence.
    """
    possibilities = 3 * reduce(mul, [
        number_of_codons[amino_acid] for amino_acid in aa_sequence
        ], 1)
    
    return(possibilities)


# In[8]:


test_possibilities = calculate_mrna_possibilities(test_sequence)

print(test_possibilities)


# So far so good. Ready to test on the real sequence!

# In[12]:


sequence = read_protein_sequence("data/rosalind_mrna.txt")

possibilities = calculate_mrna_possibilities(sequence)

answer = possibilities % 1e6

print(answer)


# Right, an error. Isn't that just what we were warned for?
# Maybe we need the modulo earlier!

# In[13]:


def calculate_mrna_possibilities_modulo(aa_sequence):
    """
Given a protein sequence (amino acids),
calculate the number of mRNA (nucleotide) sequences
that can translate into that protein sequence.
    """
    possibilities = 3 * reduce(mul, [
        number_of_codons[amino_acid] for amino_acid in aa_sequence
        ], 1) % 1e6
    
    return(possibilities)


# In[14]:


sequence = read_protein_sequence("data/rosalind_mrna.txt")

possibilities = calculate_mrna_possibilities_modulo(sequence)

print(possibilities)


# So this is not yet a solution... Let me think how to handle this then...
# 
# Right! I think I need to spend some more attention to **modular arithmetic** and **congruent modulo**.  
# I should be able to use those to solve this problem.

# Or else try this first:

# In[15]:


def calculate_mrna_possibilities_modulo2(aa_sequence):
    """
Given a protein sequence (amino acids),
calculate the number of mRNA (nucleotide) sequences
that can translate into that protein sequence.
    """
    possibilities = 3 * reduce(mul, [
        number_of_codons[amino_acid] for amino_acid in aa_sequence
        ], 1) % int(1e6)
    
    return(possibilities)


# In[17]:


print(calculate_mrna_possibilities_modulo2(sequence))


# So when the big `int` does not have to be converted to a `float`, which was necessary in this case because `1e6` is a `float` and to modulo those numbers, they need to be the same type. So that big product had to be converted. But `1e6` can be an `int` just as well, so if that is converted to an `int`, the calculation can be done and I get an answer.
# 
# This answer seems a lot easier than figuring out how to incorporate modular arithmetic. Perhaps it feels a bit like cheating, but the main goal was to make a working script, right? Not necessarily to do it by the one path that may have been intended by whoever made this exercise.
# 
# Let's go again!

# In[18]:


new_sequence = read_protein_sequence("data/rosalind_mrna2.txt")

new_possibilities = calculate_mrna_possibilities_modulo2(new_sequence)

print(new_possibilities)


# ## Success!!
# 
# That worked.
