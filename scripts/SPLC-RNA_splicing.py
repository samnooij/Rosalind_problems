#!/usr/bin/env python
# coding: utf-8

# # RNA splicing, from [Rosalind.info](https://www.rosalind.info)
# 
# (Specific exercise can be found at: rosalind.info/problems/splc/)
# 
# ## My interpretation
# 
# 1. The exercise is about splicing of mRNA molecules, and translation into amino acids
# 
# 2. The test provides a number of DNA sequences, the first of which is the sequence to be transcribed and translated. The sequences after that are introns, which are supposed to be spliced out.
# 
# Therefore, I propose to develop a script that:
# 
# 1. Reads sequences from a fasta file,
# 
# 2. Saves the first sequence as 'gene',
# 
# 3. Any other sequences are stored as introns.
# 
# 4. Then, the introns need to be found and spliced out of the 'gene' sequence
# 
# 5. and the remaining sequence is to be translated into amino acids!
# 
# (For step 5, I could insert a hand-made translation table, but it makes more sense to use Biopython's Seq objects to do the trick.)

# In[1]:


from Bio.Seq import Seq #for working with Seq objects and translating sequences
from Bio import SeqIO #for reading fasta files
from pathlib import Path #for working with files in general


# In[2]:


def extract_gene_and_introns_from_fasta(input_file):
    """
Given an input file of fasta format,
read it and save the first sequence as gene of interest,
and any other sequences as introns.
    """
    sequence_list = []
    
    for seq_record in SeqIO.parse(input_file, "fasta"):
        sequence = seq_record.seq
        #I will only need sequences for this exercise
        #Otherwise, other features may be extracted from the
        # seq_record here.
    
        sequence_list.append(sequence)

    return sequence_list


# In[3]:


test_file = Path("data/Example_RNA_splicing.fasta")

test_sequences = extract_gene_and_introns_from_fasta(test_file)

print(test_sequences)


# So far so good. Reading the sequences and saving them in a list works.
# 
# There is only one thing not entirely clear from the exercise: what if an intron sequence occurs multiple times in the gene? Should each copy be spliced? Or only the first?
# 
# Let's for now assume any number of matching introns should be spliced out.

# In[10]:


def splice_introns(gene, introns):
    """
Given a gene sequence and intron sequences,
extract the intron sequences from the gene and return
the concatenated exons.
    """
    gene = str(gene) #Using Python's string.replace() method makes splicing easier.
    
    for intron in introns:
        gene = gene.replace(str(intron), "")
        
    return Seq(gene)


# In[11]:


test_spliced_sequence = splice_introns(test_sequences[0],
                                       test_sequences[1:])

print(test_spliced_sequence)


# In[12]:


print(test_spliced_sequence.translate())


# The example on the website provides this protein sequence as solution:
# 
# `MVYIADKQHVASREAYGHMFKVCA`
# 
# So my answer is correct, but I should remove the stop codon (asterisk) from the final answer.

# In[14]:


print(test_spliced_sequence.translate().strip("*"))


# Right, that should do the trick!
# 
# Now let's jump right into the test:

# In[16]:


sequence_file = Path("data/rosalind_splc.txt")

sequence_list = extract_gene_and_introns_from_fasta(sequence_file)

spliced_sequence = splice_introns(sequence_list[0],
                                 sequence_list[1:])

translated_exons = spliced_sequence.translate().strip("*")

print(translated_exons)


# ## Success!!
# 
# I solved the exercise. That went quick!
