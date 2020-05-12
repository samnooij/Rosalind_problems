#!/usr/bin/env python
# coding: utf-8

# # Finding a shared motif, from [rosalind.info](http://rosalind.info/)
# 
# (text copied from http://rosalind.info/problems/lcsm/)
# 
# > **Problem**
# >
# > A common substring of a collection of strings is a substring of every member of the collection. We say that a common substring is a longest common substring if there does not exist a longer common substring. For example, "CG" is a common substring of "A**CG**TACGT" and "AAC**CG**TATA", but it is not as long as possible; in this case, "CGTA" is a longest common substring of "A**CGTA**CGT" and "AAC**CGTA**TA".
# >
# > Note that the longest common substring is not necessarily unique; for a simple example, "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".
# >
# > **Given**: A collection of k (k≤100) DNA strings of length at most 1 kbp each in FASTA format.
# >
# > **Return**: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)
# >
# > **Sample Dataset**
# >
# > ```
# > >Rosalind_1
# > GATTACA
# > >Rosalind_2
# > TAGACCA
# > >Rosalind_3
# > ATACA
# > ```
# >
# > **Sample Output**
# >
# > ```
# > AC
# > ```

# ## My interpretation/reasoning
# 
# 1. I am going to look for common substrings of DNA sequences in fasta format.  
#     From the example, it seems that 'common' means 'shared by _all_ sequences'!
# 
# 2. Of all the common substrings that exist, I only want to find the longest.
# 
# 3. If multiple substrings are the longest, I may return a random one. (Any correct answer is good.)
# 
# 4. There will be no more than 100 DNA sequences, with a length up to 1,000 bp.
# 
# Now practically, this means the code is to:
# 
#   - Open and read a fasta text file
#   - Take the first sequence
#   - In the second sequence, look for overlaps (>1) with the first and save all of them
#   - For each following sequence, compare to all overlapping sequences and overwrite when necessary (when the new overlap is shorter)
#   - In the end, find the longest remaining overlapping sequence. If multiple exist, pick the first.
#   
# The most difficult part will be in finding overaps/common substrings. 
# It looks like the `regex` library has nice functionality for that.
# (Found thanks to David C at https://stackoverflow.com/a/18966891).

# In[56]:


from Bio import SeqIO
import regex

def read_fasta_file(fasta_file):
    """
Takes a fasta file as input to return a dictionary with sequence IDs as keys and sequences as values.
    """
    sequence_dict = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence_dict[record.id] = record.seq
        
    return(sequence_dict)

def find_initial_overlaps(seq1, seq2):
    """
List all overlaps between the first two sequences.
    """
    substring_lengths = range(2, len(seq1))[::-1]
    #List all possible substring lengths, 
    # from the length of sequence 1 minus 1 down to 2.
    
    common_substrings = []
    #Save all common substrings in a list

    for length in substring_lengths:
        #For each possible length of substrings...
        substring_list = regex.findall(r"\w{%i}" % length, seq1, overlapped=True)
        #... list all possible substrings.
        
        for substring in substring_list:
            #For each of the substrings
            if substring in seq2 and substring not in common_substrings:
                #if it is in sequence 2 and not yet in the list
                common_substrings.append(substring)
                #add it to the list.
            else:
                #Either the substring is not shared between the sequences,
                # or it has been added to the list already.
                pass
    
    #Return all substrings common to sequences 1 and 2:
    return(common_substrings)

def compare_overlaps(motif_list, seq):
    """
Compare the current sequence (from 3 to k) to the list of overlaps (motifs) that have been found to far.
    """
    common_substrings = []
    #Keep track of the substrings that are common between the
    # motif list and the current sequence.
    
    for motif in motif_list:
        #For each motif previously identified
        if motif in seq:
            #Check if it exists in the current sequence
            common_substrings.append(motif)
            #and save it if it is.
        else:
            #Otherwise, it is no longer common between all sequences,
            # so pass. Do not track this motif any longer.
            pass
    
    #Return only substrings/motifs that are shared with the current sequence.
    return(common_substrings)

def pick_longest(motif_list):
    """
Given a list of motifs take the longest motif. 
Takes the first occurrence if multiple tie for longest, i.e. in alphabetical order.
    """
    motif_list = sorted(motif_list) #sort alphabetically first
    longest = sorted(motif_list, key = len, reverse = True)[0]
    return(longest)

def find_a_shared_motif(input_file):
    """
The complete program as one function:
 1. read the fasta file
 2. find overlaps between sequences 1 and 2 (create initial motif list)
 3. compare motif list with all other sequences, keeping only those they
    have in common
 4. return the longest motif common to all sequences
    """
    sequence_dict = read_fasta_file(input_file)
    
    initial_motifs = find_initial_overlaps(
                       str(list(sequence_dict.values())[0]),
                       str(list(sequence_dict.values())[1]))
    
    motif_list = initial_motifs
    
    for other_sequence in list(sequence_dict.values())[2:]:
        motif_list = compare_overlaps(motif_list, str(other_sequence))

    return(pick_longest(motif_list))


# In[45]:


print(find_initial_overlaps("ABCDEFG", "DEFGHI"))


# In[34]:


regex.findall(r"\w{3}", "ABCD", overlapped=True)


# In[57]:


find_a_shared_motif("data/Example_finding_a_shared_motif.txt")


# I think I have a pretty decent algorithm now. Let's try the 'real' dataset and see how that goes:

# In[58]:


find_a_shared_motif("data/rosalind_lcsm.txt")


# ## Success!!
# 
# That worked! Now let's save this notebook and commit to git.
