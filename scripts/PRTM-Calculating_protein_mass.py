#!/usr/bin/env python
# coding: utf-8

# # Calculating protein mass, from [Rosalind.info](https://www.rosalind.info)
# 
# (Specific exercise can be found at: http://rosalind.info/problems/prtm/)
# 
# ## My personal interpretation
# 
# 1. The exercise is about calculating the molecular weight of a protein
# 
# 2. The protein is represented as an amino acid sequence (a string of letters)
# 
# 3. Molecular weights per amino acid are given in a table of monoisotopic masses
# 
# 4. The practical side of the exercise comes down to reading the table with masses and then translating the letters from a given sequence into numbers using the table and adding the numbers up.
# 
# I think I can do this in three functions:
# 
# 1. Read the monoisotopic mass table and convert to a dictionary
# 
# 2. Read the text file with the amino acid sequence
# 
# 3. Take the amino acid sequence and mass table to calculate the mass

# In[5]:


def read_monoisotopic_mass_table(input_file):
    """
Given a tab-separatedd input file with amino acids (as capital letters)
in the first column, and molecular weights (as floating point numbers)
in the second column - create a dictionary with the amino acids as keys
and their respective weights as values.
    """
    mass_dict = {}
    
    with open(input_file, "r") as read_file:
        for line in read_file:
            elements = line.split()
            amino_acid = str(elements[0])
            weight = float(elements[1])
            
            mass_dict[amino_acid] = weight
            
    return mass_dict


# In[6]:


mass_dict = read_monoisotopic_mass_table("data/monoisotopic_mass_table.tsv")

print(mass_dict)


# So far so good, now make the second function:

# In[13]:


def read_amino_acid_sequence(input_file):
    """
Read a text file with an amino acid sequence and
return the sequence as string.
    """
    with open(input_file, "r") as read_file:
        for line in read_file:
            amino_acids = str(line.strip())
            #Note: the .strip() is required to remove the
            # newline, which otherwise would be interpreted
            # as amino acid!
            
    return amino_acids


# In[14]:


example_protein = read_amino_acid_sequence("data/Example_calculating_protein_mass.txt")

print(example_protein)


# Now that works as well, time to make the final function: the one that converts the amino acid sequence to its weight.

# In[15]:


def calculate_protein_weight(protein, mass_table):
    """
Given a protein sequence as string and a mass table as dictionary
(with amino acids as keys and their respective weights as values),
calculate the molecular weight of the protein by summing up the
weight of each amino acid in the protein.
    """
    total_weight = 0
    
    for amino_acid in protein:
        weight = mass_table[amino_acid]
        total_weight += weight
        
    return total_weight


# In[16]:


calculate_protein_weight(example_protein, mass_dict)


# Now this answer looks good, except the rounding of the decimals is slightly different from the example on rosalind.info... Perhaps I should just round the answer to 3 decimals?

# In[17]:


round(calculate_protein_weight(example_protein, mass_dict), 3)


# Perfect! Now let me just overwrite the function to incorporate the rounding:

# In[18]:


def calculate_protein_weight(protein, mass_table):
    """
Given a protein sequence as string and a mass table as dictionary
(with amino acids as keys and their respective weights as values),
calculate the molecular weight of the protein by summing up the
weight of each amino acid in the protein.
    """
    total_weight = 0
    
    for amino_acid in protein:
        weight = mass_table[amino_acid]
        total_weight += weight
        
    return round(total_weight, 3)


# And let's give the actual exercise a shot with this!

# In[19]:


test_protein = read_amino_acid_sequence("data/rosalind_prtm.txt")

molecular_weight = calculate_protein_weight(test_protein, mass_dict)

print(molecular_weight)

