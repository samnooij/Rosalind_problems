#!/usr/bin/env python
# coding: utf-8

# # Finding a protein motif, from [Rosalind.info](https://www.rosalind.info)
# 
# (Text copied from http://rosalind.info/problems/mprt/)
# 
# 
# <div class="problem-statement problem-statement-bordered" problem="241">
#     <blockquote>
# <h2 id="motif-implies-function">Motif Implies Function</h2>
# <div class="thumb"><a figure="Figure 1" href="http://rosalind.info/media/problems/mprt/cyclophilines.png" lightbox-title="The human cyclophilin family, as represented by the structures of the isomerase domains of some of its members." rel="lightbox[figures]"><img src="/media/problems/mprt/cyclophilines.thumb.png" /></a><div class="caption"><strong>Figure 1</strong><span>. </span><span>The human cyclophilin family, as represented by the structures of the isomerase domains of some of its members.</span></div>
# </div>
# <p>As mentioned in <a data-content="Solved by 15945 (correct ratio 67.6%)." data-trigger="hover" href="/problems/prot/" rel="popover" title="“Translating RNA into Protein”">“Translating RNA into Protein”</a>, <a class="term" href="/glossary/protein/" id="term-206" rel="tooltip" title="The functional unit of the cell.">proteins</a> perform every practical function in the <a class="term" href="/glossary/cell/" id="term-257" rel="tooltip" title="
# The &quot;building block of life,&quot; making up all living things on Earth.">cell</a>.
# A structural and functional unit of the protein is a <a class="term new" href="/glossary/protein-domain/" id="term-562" rel="tooltip" title="New term: 
# A structural and functional unit of the protein.">protein domain</a>: in terms of the protein's
# <a class="term" href="/glossary/protein-primary-structure/" id="term-570" rel="tooltip" title="
# The order of amino acids on a protein.">primary structure</a>, the domain is an interval of amino acids that can evolve and
# function independently.</p>
# <p>Each domain usually corresponds to a single function of the protein (e.g., binding the protein to <a class="term" href="/glossary/dna/" id="term-545" rel="tooltip" title="
# The molecule encoding heredity and underlying the cellular processes of all life forms.">DNA</a>, creating
# or breaking specific chemical bonds, etc.). Some proteins, such as myoglobin and the Cytochrome complex,
# have only one domain, but many proteins are multifunctional and therefore possess several domains.
# It is even possible to artificially fuse different domains into a protein molecule with definite properties,
# creating a <a class="term new" href="/glossary/chimeric-protein/" id="term-583" rel="tooltip" title="New term: 
# A protein artificially constructed from several known domains.">chimeric protein</a>.</p>
# <p>Just like species, proteins can evolve, forming <a class="term" href="/glossary/homologous/" id="term-324" rel="tooltip" title="
# Descending from the same ancestor.">homologous</a> groups called <a class="term new" href="/glossary/protein-family/" id="term-586" rel="tooltip" title="New term: 
# A group of homologous proteins.">protein families</a>.
# Proteins from one family usually have the same set of domains, performing similar functions;
# see <a href="/media/problems/mprt/cyclophilines.png" lightbox-title="The human cyclophilin family, as represented by the structures of the isomerase domains of some of its members." rel="lightbox[figures]" title="Click to view">Figure 1</a>.</p>
# <p>A component of a domain essential for its function is called a <a class="term" href="/glossary/motif/" id="term-241" rel="tooltip" title="
# A nucleotide or amino acid pattern of biological significance.">motif</a>, a term that in general
# has the same meaning as it does in <a class="term" href="/glossary/nucleic-acid/" id="term-200" rel="tooltip" title="
# A polymer of nucleotides, constituting either RNA or DNA.">nucleic acids</a>, although many other terms are also used
# (blocks, signatures, fingerprints, etc.) Usually protein motifs are evolutionarily conservative,
# meaning that they appear without much change in different species.</p>
# <p>Proteins are identified in different labs around the world and gathered into freely accessible databases.
# A central repository for protein data is <a href="http://www.uniprot.org/" target="_blank">UniProt</a>, which provides
# detailed protein annotation, including function description, domain structure, and post-translational modifications.
# UniProt also supports protein similarity search, taxonomy analysis, and literature citations.</p>
# </blockquote>
# <h2 id="problem">Problem</h2>
# <p>To allow for the presence of its varying forms, a protein motif is represented by a shorthand as follows:
# [XY] means "either X or Y" and {X} means "any amino acid except X."  For example, the N-glycosylation motif
# is written as N{P}[ST]{P}.</p>
# <p>You can see the complete description and features of a particular protein by its access ID
# "uniprot_id" in the UniProt database, by inserting the ID number into</p>
# <div class="codehilite"><pre>http://www.uniprot.org/uniprot/uniprot_id
# </pre></div>
# 
# <p>Alternatively, you can obtain a protein sequence in <a class="term" href="/glossary/fasta-format/" id="term-759" rel="tooltip" title="
# A text format used for naming genetic strings in databases.">FASTA format</a> by following</p>
# <div class="codehilite"><pre>http://www.uniprot.org/uniprot/uniprot_id.fasta
# </pre></div>
# 
# 
# <p>For example, the data for protein B5ZC00 can be found at <a href="http://www.uniprot.org/uniprot/B5ZC00" target="_blank"><a href="http://www.uniprot.org/uniprot/B5ZC00" rel="nofollow" target="_blank"><a href="http://www.uniprot.org/uniprot/B5ZC00" rel="nofollow" target="_blank">http://www.uniprot.org/uniprot/B5ZC00</a></a></a>.</p>
# <p><span class="given-return">Given:</span> At most 15 UniProt Protein Database access IDs.</p>
# <p><span class="given-return">Return:</span> For each protein possessing the N-glycosylation motif, output its given access ID followed
# by a list of <a class="term" href="/glossary/location/" id="term-382" rel="tooltip" title="
# The position in a string where a substring begins.">locations</a> in the protein string where the motif can be found.</p>
# <h2 id="sample-dataset">Sample Dataset</h2>
# <div class="codehilite"><pre>A2Z669
# B5ZC00
# P07204_TRBM_HUMAN
# P20840_SAG1_YEAST
# </pre></div>
# 
# 
# <h2 id="sample-output">Sample Output</h2>
# <div class="codehilite"><pre>B5ZC00
# 85 118 142 306 395
# P07204_TRBM_HUMAN
# 47 115 116 382 409
# P20840_SAG1_YEAST
# 79 109 135 248 306 348 364 402 485 501 614
# </pre></div>
# 
# 
# <blockquote>
# <h2 id="note">Note</h2>
# <p>Some entries in UniProt have one primary (citable) accession number and some secondary numbers, appearing due to
# merging or demerging entries. In this problem, you may be given any type of ID.
# If you type the secondary ID into the UniProt query, then you will be automatically
# redirected to the page containing the primary ID.
# You can find more information about UniProt IDs <a href="http://www.uniprot.org/manual/accession_numbers" target="_blank">here</a>.</p>
# </blockquote>
# <div class="clearfix"></div>
# </div>

# ## My interpretation/reasoning
# 
# 1. In this exercise, I will combine two things:
#     1. Downloading protein sequences from a public database
#     2. Finding motifs in the sequences
#     
# 2. The motif of interest is "N, not-P, S-or-T, not-P" so a 4-amino acid sequence.
# 
# 3. The results should be the protein ID of the protein that contains the motif, followed by a newline and all positions of the motif, separated by spaces.
# 
# So practically, I want to make a script that:
#   - Opens and reads a text file with IDs
#   - For each ID, lookup the amino acid sequence
#   - Find any position that holds the motif N{P}[ST]{P}
#       - if there are none: pass
#       - if the motif is found: return the ID and the positions
#       
# Sounds pretty straightforward. Let's see how to get that into code.

# In[3]:


def read_ids_from_file(input_file):
    """
Given a text file with one ID per line, 
return a list of IDs.
    """
    id_list = []
    
    with open(input_file, 'r') as read_file:
        for line in read_file:
            id_list.append(line.strip())
            
    return(id_list)


# In[3]:


print(read_ids_from_file("data/Example_finding_a_protein_motif.txt"))


# Then I want to download the fasta files belonging to these IDs, for which I can use the urllib module as suggested here: https://stackoverflow.com/questions/1393324/in-python-given-a-url-to-a-text-file-what-is-the-simplest-way-to-read-the-cont

# In[6]:


import urllib.request


# In[5]:


def fetch_uniprot_fasta(accession_id):
    """
Given an accession ID from UniProt,
download and return the corresponding fasta record
as separate variables for ID and sequence.
    """
    base_url = "http://www.uniprot.org/uniprot/"
    
    fasta_url = base_url + accession_id  + ".fasta"
    
    sequence = ""
    
    for line in urllib.request.urlopen(fasta_url):
        text = line.decode("utf-8").strip()
        
        if text.startswith(">"):
            fasta_id = text
        else:
            sequence += text
        
    return(fasta_id, sequence)


# In[29]:


fetch_uniprot_fasta("A2Z669")


# So this way you can just read and print a fasta file from the web. Can we also immediately parse it as fasta record?

# In[30]:


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# In[31]:


def create_seq_record(name, sequence):
    """
Given a name and a sequence,
create a Biopython SeqRecord object.
    """
    sequence = Seq(sequence)
    seq_record = SeqRecord(sequence)
    seq_record.id = name
    
    return(seq_record)


# In[34]:


(name, sequence) = fetch_uniprot_fasta("A2Z669")

print(create_seq_record(name, sequence))


# This is not strictly necessary, but I thought it was nice to have. 
# Let's just continue with looking for motifs in the sequences.

# In[8]:


import regex


# In[72]:


def search_motif_positions(sequence, motif):
    """
Search a sequence for a given motif (as regex),
and return position as list. Return False if the motif
is absent from the sequence.
    """
    positions = [m.start() + 1 for m in regex.finditer(
        r"%s" % motif, sequence)]
    #Method thanks to moinudin: https://stackoverflow.com/a/4664889
    
    if len(positions) > 0:
        return(positions)
    else:
        return(False)


# The regex of the motif should be: "N[^P][S|T][^P]"

# In[73]:


print(search_motif_positions(sequence, "N[^P][S|T][^P]"))


# In[64]:


regex.findall(r"N[^P][S|T][^P]", "MKNKFKTQEELVNHLKTVGFVFANSEIYNGLANAWDYGPLGVLLKNNLKNLWWKEFVTKQKDVVGLDSAIILNPLVWKASGHLDNFSDPLIDCKNCKARYRADKLIESFDENIHIAENSSNEEFAKVLNDYEISCPTCKQFNWTEIRHFNLMFKTYQGVIEDAKNVVYLRPETAQGIFVNFKNVQRSMRLHLPFGIAQIGKSFRNEITPGNFIFRTREFEQMEIEFFLKEESAYDIFDKYLNQIENWLVSACGLSLNNLRKHEHPKEELSHYSKKTIDFEYNFLHGFSELYGIAYRTNYDLSVHMNLSKKDLTYFDEQTKEKYVPHVIEPSVGVERLLYAILTEATFIEKLENDDERILMDLKYDLAPYKIAVMPLVNKLKDKAEEIYGKILDLNISATFDNSGSIGKRYRRQDAIGTIYCLTIDFDSLDDQQDPSFTIRERNSMAQKRIKLSELPLYLNQKAHEDFQRQCQK")


# In[ ]:





# In[82]:


def find_a_protein_motif(input_file):
    """
The complete program:
 - read a file with accession IDs
 - parse the associated amino acid sequences
 - look for the motif
 - if the motif exists:
   - print accession ID and (start) positions
    """
    motif = "N[^P][S|T][^P]"
    
    accession_ids = read_ids_from_file(input_file)
    
    for accession_id in accession_ids:
        (fasta_id, sequence) = fetch_uniprot_fasta(accession_id)
        
        if search_motif_positions(sequence, motif):
            print("%s\n%s" % (accession_id, " ".join(map(str, search_motif_positions(sequence, motif)))))
        else:
            pass
        
    return(None)


# In[83]:


find_a_protein_motif("data/Example_finding_a_protein_motif.txt")


# Now that seems to be working. Let's see if it works on the real deal:

# In[84]:


find_a_protein_motif("data/rosalind_mprt.txt")


# Apparently this was the wrong answer, unfortunately. Where could the error be?
# 
# What if I just try it again? Maybe it was just this set for whatever reason.

# In[85]:


find_a_protein_motif("data/rosalind_mprt2.txt")


# No, that didn't do it. Then I have to find something else...
# 
# Let's make a test version to debug the problem. A function that shows all IDs, sequences, motifs, positions to help check if it all looks right.

# In[1]:


def find_a_protein_motif_debug(input_file):
    """
The complete program:
 - read a file with accession IDs
 - parse the associated amino acid sequences
 - look for the motif
 - if the motif exists:
   - print accession ID and (start) positions
    """
    motif = "N[^P][S|T][^P]"
    
    accession_ids = read_ids_from_file(input_file)
    
    for accession_id in accession_ids:
        print("ID: %s" % accession_id)
        (fasta_id, sequence) = fetch_uniprot_fasta(accession_id)
        
        print("Fasta ID: %s\nSequence: %s" % (fasta_id, sequence))
        
        search_motif_positions_debug(sequence, motif)
        
    return(None)

def search_motif_positions_debug(sequence, motif):
    """
Search a sequence for a given motif (as regex),
and return position as list. Return False if the motif
is absent from the sequence.
    """
    
    print("Found motifs: %s" % regex.findall(r"%s" % motif, sequence))
    
    positions = [m.start() + 1 for m in regex.finditer(
        r"%s" % motif, sequence)]
    #Method thanks to moinudin: https://stackoverflow.com/a/4664889
    
    print("Motif positions: %s" % positions)
    
    if len(positions) > 0:
        return(positions)
    else:
        return(False)


# In[9]:


find_a_protein_motif_debug("data/rosalind_mprt.txt")


# When looking back at the example sequences, I found the problem! Consider these example answers:
# 
# **By Rosalind**
# ```
# B5ZC00
# 85 118 142 306 395
# P07204_TRBM_HUMAN
# 47 115 116 382 409
# P20840_SAG1_YEAST
# 79 109 135 248 306 348 364 402 485 501 614
# ```
# 
# **By my function**
# ```
# B5ZC00
# 85 118 142 306 395
# P07204_TRBM_HUMAN
# 47 115 382 409
# P20840_SAG1_YEAST
# 79 109 135 248 306 348 364 402 485 501 614
# ```
# 
# And look carefully at the list of positions given for "P07204_TRBM_HUMAN"...
# 
# ... see it yet?
# 
# The answer is...
# 
# `116` is missing from my function's output! It does not list overlapping matches. 
# So I need to find a regex function that matches all overlapping matches too!
# 
# Perhaps take another look at: https://stackoverflow.com/a/18966891?

# In[21]:


#Let's take this sequence specifically and see how to fix my problem there!

test_id = ">sp|P07204|TRBM_HUMAN Thrombomodulin OS=Homo sapiens OX=9606 GN=THBD PE=1 SV=2"

test_sequence = "MLGVLVLGALALAGLGFPAPAEPQPGGSQCVEHDCFALYPGPATFLNASQICDGLRGHLMTVRSSVAADVISLLLNGDGGVGRRRLWIGLQLPPGCGDPKRLGPLRGFQWVTGDNNTSYSRWARLDLNGAPLCGPLCVAVSAAEATVPSEPIWEEQQCEVKADGFLCEFHFPATCRPLAVEPGAAAAAVSITYGTPFAARGADFQALPVGSSAAVAPLGLQLMCTAPPGAVQGHWAREAPGAWDCSVENGGCEHACNAIPGAPRCQCPAGAALQADGRSCTASATQSCNDLCEHFCVPNPDQPGSYSCMCETGYRLAADQHRCEDVDDCILEPSPCPQRCVNTQGGFECHCYPNYDLVDGECVEPVDPCFRANCEYQCQPLNQTSYLCVCAEGFAPIPHEPHRCQMFCNQTACPADCDPNTQASCECPEGYILDDGFICTDIDECENGGFCSGVCHNLPGTFECICGPDSALARHIGTDCDSGKVDGGDSGSGEPPPSPTPGSTLTPPAVGLVHSGLLIGISIASLCLVVALLALLCHLRKKQGAARAKMEYKCAAPSKEVVLQHVRTERTPQRL"

test_matches = regex.findall(r"N[^P][S|T][^P]", test_sequence, overlapped=True)

print([match for match in test_matches])

test_iters = regex.finditer(r"N[^P][S|T][^P]", test_sequence, overlapped=True)

print([match.start() + 1 for match in test_iters])


# It looks like I was missing only the `overlapped=True` part. If I add that, I should be good.

# In[22]:


def search_motif_positions2(sequence, motif):
    """
Search a sequence for a given motif (as regex),
and return position as list. Return False if the motif
is absent from the sequence.
    """
    positions = [m.start() + 1 for m in regex.finditer(
        r"%s" % motif, sequence, overlapped = True)]
    #Method thanks to moinudin: https://stackoverflow.com/a/4664889
    
    if len(positions) > 0:
        return(positions)
    else:
        return(False)


# In[23]:


def find_a_protein_motif2(input_file):
    """
The complete program:
 - read a file with accession IDs
 - parse the associated amino acid sequences
 - look for the motif
 - if the motif exists:
   - print accession ID and (start) positions
    """
    motif = "N[^P][S|T][^P]"
    
    accession_ids = read_ids_from_file(input_file)
    
    for accession_id in accession_ids:
        (fasta_id, sequence) = fetch_uniprot_fasta(accession_id)
        
        if search_motif_positions2(sequence, motif):
            print("%s\n%s" % (accession_id, " ".join(map(str, search_motif_positions2(sequence, motif)))))
        else:
            pass
        
    return(None)


# So now these should give different results than before if I use the same input files. Let's see how that goes.

# In[24]:


find_a_protein_motif2("data/Example_finding_a_protein_motif.txt")


# So far so good...

# In[25]:


find_a_protein_motif2("data/rosalind_mprt.txt")


# Indeed there is a minor difference here! Look at the positions `168` and `169` for `P01047_KNL2_BOVIN`.

# In[26]:


find_a_protein_motif2("data/rosalind_mprt2.txt")


# So this one should differ in the positions for `P02974_PMM1_NEIGO`: `68` has been added compared to the previous attempt.
# 
# Now I feel I am ready for the challenge again. Let's do it!

# In[27]:


find_a_protein_motif2("data/rosalind_mprt3.txt")


# ## Success!!
# 
# I did it! Indeed the overlapping motifs were what I missed first and what has been fixed in this second version.
