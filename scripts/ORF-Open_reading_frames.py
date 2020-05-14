#!/usr/bin/env python
# coding: utf-8

# # Open reading frames, from [Rosalind.info](https://www.rosalind.info)
# 
# (Text copied from http://rosalind.info/problems/orf/)
# 
# <div class="problem-statement problem-statement-bordered" problem="122">
# <h2 id="problem">Problem</h2>
# <p>Either strand of a DNA double helix can serve as the <a class="term" href="/glossary/coding-strand/" id="term-322" rel="tooltip" title="
# The strand of a double-stranded DNA molecule that is copied or transcribed into RNA.">coding strand</a> for RNA transcription.
# Hence, a given DNA string implies six total <a class="term new" href="/glossary/reading-frame/" id="term-267" rel="tooltip" title="New term: 
# One of three possible ways to read a given strand of DNA, depending upon the starting position.">reading frames</a>, or ways in which the same region of DNA can be translated into amino acids:
# three reading frames result from reading the string itself, whereas three more result
# from reading its <a class="term" href="/glossary/reverse-complement/" id="term-252" rel="tooltip" title="
# The DNA string formed by reversing and complementing each symbol.">reverse complement</a>.</p>
# <p>An <a class="term new" href="/glossary/open-reading-frame/" id="term-268" rel="tooltip" title="New term: 
# A sequence in DNA or RNA potentially able to encode the protein.">open reading frame</a> (ORF) is one which starts from the <a class="term" data-math="true" href="/glossary/start-codon/" id="term-260" rel="tooltip" title="
# The RNA codon $\textrm{AUG}$, which codes for the amino acid methionine and indicates the beginning
# of translation into protein.">start codon</a> and ends by <a class="term" href="/glossary/stop-codon/" id="term-261" rel="tooltip" title="
# One of three possible RNA codons that indicate the termination of protein translation.">stop codon</a>, without
# any other <a class="term" href="/glossary/stop-codon/" id="term-261" rel="tooltip" title="
# One of three possible RNA codons that indicate the termination of protein translation.">stop codons</a> in between. Thus, a candidate protein string is derived by translating an open reading
# frame into amino acids until a stop codon is reached.</p>
# <p><span class="given-return">Given:</span> A <a class="term" href="/glossary/dna-string/" id="term-349" rel="tooltip" title="
# A string constructed from the alphabet {A, C, G, T}.">DNA string</a> <mathjax>$s$</mathjax> of length at most 1 <a class="term" href="/glossary/kbp/" id="term-394" rel="tooltip" title="
# 1 kbp = 1000 base pairs">kbp</a> in <a class="term" href="/glossary/fasta-format/" id="term-759" rel="tooltip" title="
# A text format used for naming genetic strings in databases.">FASTA format</a>.</p>
# <p><span class="given-return">Return:</span> Every distinct candidate protein string that can be translated from ORFs of <mathjax>$s$</mathjax>.
# Strings can be returned in any order.</p>
# <h2 id="sample-dataset">Sample Dataset</h2>
# <div class="codehilite"><pre>&gt;Rosalind_99
# AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG
# </pre></div>
# 
# 
# <h2 id="sample-output">Sample Output</h2>
# <div class="codehilite"><pre>MLLGSFRLIPKETLIQVAGSSPCNLS
# M
# MGMTPRLGLESLLE
# MTPRLGLESLLE
# </pre></div>

# ## My interpretation/reasoning
# 
# 1. This problem seems straightforward: return all possible six frame translations for a given DNA sequence
# 
# 2. Several tools exist that predict and translate open reading frames, such as [this useful webtool by ExPASy](https://web.expasy.org/translate/), or [Prodigal](https://github.com/hyattpd/Prodigal) also seems to use something like six frame translations. 
# 
# 3. However, for this particular exercise, these tools do not seem so handy: the  output needs to be formatted as above: One translation per line, each possible translation is requested, and nothing else.
# 
# So...
# 
# To put this into a script, I want to have:
#  - A function to read the sequence from a fasta file
#  - Translations for the DNA sequences in all 6 possible open reading frames, _not_ stopping at any one particular stop codon. That is, stop the current translation, but keep looking in the rest of the sequence.
#  - This bit of Biopython's documentation can help me do this: https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc295
#  
# In particular, this piece of code may be of use:
# 
# ```python
# >>> from Bio import SeqIO
# >>> record = SeqIO.read("NC_005816.fna", "fasta")
# >>> table = 11
# >>> min_pro_len = 100
# 
# Here is a neat trick using the Seq object’s split method to get a list of all the possible ORF translations in the six reading frames:
# 
# >>> for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
# ...     for frame in range(3):
# ...         length = 3 * ((len(record)-frame) // 3) #Multiple of three
# ...         for pro in nuc[frame:frame+length].translate(table).split("*"):
# ...             if len(pro) >= min_pro_len:
# ...                 print("%s...%s - length %i, strand %i, frame %i" \
# ...                       % (pro[:30], pro[-3:], len(pro), strand, frame))
# GCLMKKSSIVATIITILSGSANAASSQLIP...YRF - length 315, strand 1, frame 0
# KSGELRQTPPASSTLHLRLILQRSGVMMEL...NPE - length 285, strand 1, frame 1
# GLNCSFFSICNWKFIDYINRLFQIIYLCKN...YYH - length 176, strand 1, frame 1
# VKKILYIKALFLCTVIKLRRFIFSVNNMKF...DLP - length 165, strand 1, frame 1
# NQIQGVICSPDSGEFMVTFETVMEIKILHK...GVA - length 355, strand 1, frame 2
# RRKEHVSKKRRPQKRPRRRRFFHRLRPPDE...PTR - length 128, strand 1, frame 2
# TGKQNSCQMSAIWQLRQNTATKTRQNRARI...AIK - length 100, strand 1, frame 2
# QGSGYAFPHASILSGIAMSHFYFLVLHAVK...CSD - length 114, strand -1, frame 0
# IYSTSEHTGEQVMRTLDEVIASRSPESQTR...FHV - length 111, strand -1, frame 0
# WGKLQVIGLSMWMVLFSQRFDDWLNEQEDA...ESK - length 125, strand -1, frame 1
# RGIFMSDTMVVNGSGGVPAFLFSGSTLSSY...LLK - length 361, strand -1, frame 1
# WDVKTVTGVLHHPFHLTFSLCPEGATQSGR...VKR - length 111, strand -1, frame 1
# LSHTVTDFTDQMAQVGLCQCVNVFLDEVTG...KAA - length 107, strand -1, frame 2
# RALTGLSAPGIRSQTSCDRLRELRYVPVSL...PLQ - length 119, strand -1, frame 2
# ```

# In[1]:


from Bio import SeqIO


# In[3]:


def read_single_sequence_file(input_file):
    """
Read a single sequence from a file.
    """
    return(SeqIO.read(input_file, "fasta"))


# In[5]:


def translate_six_frames(seq_record):
    """
Given a Biopython sequence record with a DNA (or RNA?) sequence,
translate into amino acid (protein) sequences in six frames.
Returns translations as list of strings.
    """
    translation_list = []
    
    for strand, nuc in [(+1, seq_record.seq), (-1,  seq_record.seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(seq_record)-frame) // 3)
            for pro in nuc[frame:frame+length].translate().split("*"):
                translation_list.append(pro)
                
    return(translation_list)


# In[6]:


test_record = read_single_sequence_file("data/Example_open_reading_frames.txt")


# In[7]:


test_list = translate_six_frames(test_record)


# In[8]:


for translation in test_list: print(translation)


# This is not quite what I expected... The list is way too long and does not always start with start codons (or Ms in the translation).
# Let's debug the function a bit by using print statements.

# In[9]:


def test_translate_six_frames(seq_record):
    """
Given a Biopython sequence record with a DNA (or RNA?) sequence,
translate into amino acid (protein) sequences in six frames.
Returns translations as list of strings.
    """
    translation_list = []
    
    for strand, nuc in [(+1, seq_record.seq), (-1,  seq_record.seq.reverse_complement())]:
        print("Strand: %s\nNuc: %s" % (strand, nuc))
        for frame in range(3):
            print("Frame: %s" % frame)
            length = 3 * ((len(seq_record)-frame) // 3)
            print("Length: %s" % length)
            print("Possible translations: %s" % nuc[frame:frame+length].translate())
            for pro in nuc[frame:frame+length].translate().split("*"):
                translation_list.append(pro)
                
    return(translation_list)


# In[10]:


for translation in test_translate_six_frames(test_record): print(translation)


# I think I will need to do some searching, for example with Regular Expressions (`regex`), to filter out the open reading frames, which start with start codons and end with stop codons.

# In[11]:


import regex


# In[41]:


#Let's try a little test with the RegEx: ATG([ACGT]{3})*(TAA|TAG|TGA)
# This should start with "ATG" (start), followed by any number of triplets consisting of A, C, G and/or T, and end with one of TAA, TAG or TGA (stop codons)

test_forward = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"

test_reverse = "CTGAGATGCTACTCGGATCATTCAGGCTTATTCCAAAAGAGACTCTAATCCAAGTCGCGGGGTCATCCCCATGTAACCTGAGTTAGCTACATGGCT"

pattern = "ATG([ACGT]{3})*?(TAA|TAG|TGA)"

#Explanation:
# 1. Start with "ATG" literally (start codon)
# 2.1. Then match any number of triplets of A, C, G and/or G
# 2.2. However, since '*' is a 'greedy' match, it will look for the longest possible match
#     and I want to stop at the first stop codon (= point 3). Therefore, I added the '?'
#     to make the matching method 'lazy'. (Also see https://regex101.com/)
# 3. Stop the match at either "TAA", "TAG" or "TGA" (stop codons)

print(regex.findall(r"%s" % pattern, test_forward, overlapped=True))

print(regex.findall(r"%s" % pattern, test_reverse, overlapped=True))


# In[15]:


get_ipython().run_line_magic('pinfo2', 'regex.findall')


# In[17]:


for match in regex.finditer(r"%s" % pattern, test_forward, overlapped=True): print(match)


# In[18]:


for match in regex.finditer(r"%s" % pattern, test_reverse, overlapped=True): print(match)


# This looks a bit better. This would give 5 ORFs, which is almost consistent with the example answer, which has 4 sequences.
# 
# Let's see what happens if I translate these:

# In[20]:


from Bio.Seq import Seq


# In[43]:


print(test_forward)
print()

for match in regex.finditer(r"%s" % pattern, test_forward, overlapped=True):
    test_sequence = Seq(match.captures()[0])
    print(test_sequence)
    print(test_sequence.translate().strip("*"))


# Yes, that looks pretty nice. Now on the other strand:

# In[44]:


for match in regex.finditer(r"%s" % pattern, test_reverse, overlapped=True):
    test_sequence = Seq(match.captures()[0])
    print(test_sequence)
    print(test_sequence.translate().strip("*"))


# That also looks pretty good. Now I probably want to throw all these sequences into a set, to remove duplicates. (The example shows only one time 'M'.) So the final function should look something like...:

# In[ ]:


pattern = "ATG([ACGT]{3})*?(TAA|TAG|TGA)"

#Explanation:
# 1. Start with "ATG" literally (start codon)
# 2.1. Then match any number of triplets of A, C, G and/or G
# 2.2. However, since '*' is a 'greedy' match, it will look for the longest possible match
#     and I want to stop at the first stop codon (= point 3). Therefore, I added the '?'
#     to make the matching method 'lazy'. (Also see https://regex101.com/)
# 3. Stop the match at either "TAA", "TAG" or "TGA" (stop codons)


# In[55]:


def translate_six_frames2(seq_record):
    """
Given a Biopython sequence record with a DNA (or RNA?) sequence,
translate into amino acid (protein) sequences in six frames.
Returns translations as (deduplicated) list of strings.
    """
    translation_list = []
    
    for strand, nuc in [(+1, seq_record.seq), (-1,  seq_record.seq.reverse_complement())]:
        for match in regex.finditer(r"%s" % pattern, str(nuc), overlapped=True):
            sequence = Seq(match.captures()[0])
            translation = sequence.translate().strip("*")
            
            translation_list.append(translation)
            
    return(list(set(translation_list)))


# In[56]:


test_list2 = translate_six_frames2(test_record)


# In[57]:


for translation in test_list2: print(translation)


# Well, let's see if this is enough to pass the test.

# In[58]:


seq_record = read_single_sequence_file("data/rosalind_orf.txt")

translation_list = translate_six_frames2(seq_record)

for translation in translation_list:
    print(translation)


# ## Success!!
# 
# It worked. I passed the test.
