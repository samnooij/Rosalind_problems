#!/usr/bin/env python
# coding: utf-8

# # Locating restriction sites, from [Rosalind.info](https://www.rosalind.info)
# 
# (Specific exercise can be found at: http://rosalind.info/problems/revp/)
# 
# ## My personal interpretation of the problem
# 
# 1. The problem is about DNA palindromes (or possible restriction sites)
# 
# 2. The input will be a DNA sequence of at most 1 kbp long
# 
# 3. What the program than needs to do is read the sequence from start do end,
# 
#     - while checking for each substring of length 2-6 (the forward half of the potential 4-12 bp palindrome) if the next substring of equal length is its reverse complement
#     
# 4. The output should then be 1) the position of this site, and 2) the length of the palindrome
# 
#     - so not only keep the nucleotides, but keep the position and lengths (too)!
#     

# To tackle this problem, I want to make functions that:
# 
# 1. Read a sequence from a fasta file
# 
# 2. Find palindromes! (includes: read sequence, check reverse complements of substrings, keep positions)

# In[4]:


from Bio.Seq import Seq #for working with Seq objects and making reverse complements
from Bio import SeqIO #for reading fasta files
from pathlib import Path #for working with files in general


# In[5]:


#Set minimum and maximum lengths of palindrome halves
min_length = 2
max_length = 6


# In[6]:


def read_fasta_file(input_file):
    """
Given an input fasta file, read the fasta record (ID and sequence).
(Note: read only one sequence!)
    """
    for seq_record in SeqIO.parse(input_file, "fasta"):
        name = seq_record.id #this is not necessary in this exercise
        sequence = seq_record.seq
        break
    
    return name, sequence


# In[7]:


example_file = Path("data/Example_Locating_restriction_sites.fasta")

(example_name, example_sequence) = read_fasta_file(example_file)

print(example_sequence)


# Alright, the fasta reading function seems to be working. Now move on to the palindrome finder.  
# ...  
# After doing some quick tests with the `.reverse_complement()` method!

# In[5]:


print(example_sequence.reverse_complement())


# In[6]:


print(example_sequence[3:8].reverse_complement())


# Right, that works as I expected. Now let's try and put that into the palindrome finder function.

# In[7]:


def find_DNA_palindromes(sequence):
    """
Given a DNA sequence, find all palindromes of lengths 4-12
and report their position and length, and include the sequence
for easier debugging.

Save palindromes as a dictionary, using positions as key and
length and sequence as values (as tuple).
    """
    palindrome_dict = {}
    
    for position in range(len(sequence)):
        for length in range(2,6):
            potential_forward_half = sequence[position:position+length+1]
            potential_reverse_half = sequence[position+length+1:position+(length*2)]
            
            print("Length: %i, forward: %s, reverse: %s" % (
            length,
            potential_forward_half,
            potential_reverse_half
            )
                 )
            
            
            if potential_forward_half == potential_reverse_half.reverse_complement():
                #This is a palindrome!
                palindrome_dict[position] = (
                    length, 
                    potential_forward_half + potential_reverse_half
                )
            else:
                pass
            
        return palindrome_dict
            


# In[8]:


example_palindromes = find_DNA_palindromes(example_sequence)


# Right, so this does not seem to be working correctly yet.
# The forward sequence starts as one too long,
# the reverse sequence is too short,
# and the palindromes are returned after checking one position.
# 
# Let's see if I can fix those errors:

# In[9]:


def find_DNA_palindromes(sequence):
    """
Given a DNA sequence, find all palindromes of lengths 4-12
and report their position and length, and include the sequence
for easier debugging.

Save palindromes as a dictionary, using positions as key and
length and sequence as values (as tuple).
    """
    palindrome_dict = {}
    
    for position in range(len(sequence) + 1):
        for length in range(2,6):
            potential_forward_half = sequence[position:position+length]
            potential_reverse_half = sequence[position+length+1:position+(length*2)+1]
            
            print("Length: %i, forward: %s, reverse: %s" % (
            length,
            potential_forward_half,
            potential_reverse_half
            )
                 )
            
            
            if potential_forward_half == potential_reverse_half.reverse_complement():
                #This is a palindrome!
                print("%s and %s appear to form a palindrome!" % (
                potential_forward_half,
                potential_reverse_half
                )
                     )
                
                palindrome_dict[position] = (
                    length, 
                    potential_forward_half + potential_reverse_half
                )
            else:
                pass
            
    return palindrome_dict


# In[10]:


example_palindromes = find_DNA_palindromes(example_sequence)


# Right, I need a few more checks in here:
# 
# 1. it may be nice to print the position as well
# 
# 2. the function needs to stop when at the end of the sequence (e.g. not try and make palindromes of non-existing/empty sequences!)
# 
# Also, the reverse sequence seems to start one too late.

# In[28]:


def find_DNA_palindromes(sequence, debug=False):
    """
Given a DNA sequence, find all palindromes of lengths 4-12
and report their position and length, and include the sequence
for easier debugging.

Save palindromes as a dictionary, using positions as key and
length and sequence as values (as tuple).
    """
    #initialise empty dictionary,
    # to be filled with the palindromes
    palindrome_dict = {}
    
    #loop over the complete sequence, by positions
    for position in range(len(sequence) + 1):
        #then for each position, check for potential palindromes
        for length in range(min_length, max_length + 1):
            
            if len(sequence[position:position+length]) == length:
                #Check if the substring to check is as long as intended,
                # i.e. does not 'go over' the end of the sequence
                potential_forward_half = sequence[position:position+length]
            else:
                #If it does, break (stop checking)
                break
                
            if len(sequence[position+length:position+(length*2)]) == length:
                potential_reverse_half = sequence[position+length:position+(length*2)]
            else:
                break
            
            if debug:
                print("Position: %i, Length: %i, forward: %s, reverse: %s" % (
                    position + 1,
                    length * 2,
                    potential_forward_half,
                    potential_reverse_half
                )
                     )
            
            
            if potential_forward_half == potential_reverse_half.reverse_complement():
                #This is a palindrome!
                if debug:
                    print("%s and %s appear to form a palindrome!" % (
                    potential_forward_half,
                    potential_reverse_half
                    )
                         )
                
                palindrome_dict[position + 1] = ( #add 1 to compensate for Python's 0-based counting
                    2 * length, #multiply by 2 to give the length of the complete palindrome
                    potential_forward_half + potential_reverse_half
                )
            else:
                pass
            
    return palindrome_dict


# In[29]:


example_palindromes = find_DNA_palindromes(example_sequence, debug=True)


# In[30]:


example_palindromes = find_DNA_palindromes(example_sequence)

print(example_palindromes)


# Now that this seems to work too, let's include a final function that handles printing in the right format.

# In[31]:


def print_palindromes(palindrome_dict):
    """
Given a dictionary with palindrome positions as keys, and
lengths as first element of the value,
print the positions and lengths separated by a whitespace,
one pair per line.
    """
    for key, value in palindrome_dict.items():
        print(key, value[0])
        
    return None


# In[32]:


print_palindromes(example_palindromes)


# With that, I think I have everything covered. Let's try the actual exercise!

# In[37]:


test_file = Path("data/rosalind_revp.txt")

(test_name, test_sequence) = read_fasta_file(test_file)

palindromes = find_DNA_palindromes(test_sequence)

print_palindromes(palindromes)


# Apparently, this answer is wrong. Let's check where it may have gone wrong...

# In[38]:


print("Length of test sequence: %i" % len(test_sequence))


# In[39]:


find_mistakes = find_DNA_palindromes(test_sequence, debug = True)


# It looks alright at first glance... What may have gone wrong??
# 
# Ah, I think I get it! If a position has multiple palindromes, I can only report one because I'm using a dictionary. And dictionaries may only have each key once. So to be able to report multiple palindromes on the same position, I have to save them in a different way. E.g. in a list!

# In[27]:


def find_DNA_palindromes(sequence, debug=False):
    """
Given a DNA sequence, find all palindromes of lengths 4-12
and report their position and length, and include the sequence
for easier debugging.

Save palindromes as a list, using positions, length and
sequence as values (as tuple).
    """
    #initialise empty list to be filled with the palindromes
    palindrome_list = []
    
    #loop over the complete sequence, by positions
    for position in range(len(sequence) + 1):
        #then for each position, check for potential palindromes
        for length in range(min_length, max_length + 1):
            
            if len(sequence[position:position+length]) == length:
                #Check if the substring to check is as long as intended,
                # i.e. does not 'go over' the end of the sequence
                potential_forward_half = sequence[position:position+length]
            else:
                #If it does, break (stop checking)
                break
                
            if len(sequence[position+length:position+(length*2)]) == length:
                potential_reverse_half = sequence[position+length:position+(length*2)]
            else:
                break
            
            if debug:
                print("Position: %i, Length: %i, forward: %s, reverse: %s" % (
                    position + 1,
                    length * 2,
                    potential_forward_half,
                    potential_reverse_half
                )
                     )
            
            
            if potential_forward_half == potential_reverse_half.reverse_complement():
                #This is a palindrome!
                if debug:
                    print("%s and %s appear to form a palindrome!" % (
                    potential_forward_half,
                    potential_reverse_half
                    )
                         )
                
                palindrome_list.append((
                    position + 1, #add 1 to compensate for Python's 0-based counting
                    2 * length, #multiply by 2 to give the length of the complete palindrome
                    potential_forward_half + potential_reverse_half
                ))
            else:
                pass
            
    return palindrome_list


# In[1]:


def print_palindromes_from_list(palindrome_list):
    """
Given a list with palindrome positions, lengths and
sequences, print the positions and lengths separated by a whitespace,
one pair per line.
    """
    for palindrome in palindrome_list:
        print("%s %s" % (
            palindrome[0],
        palindrome[1]
        ))
        
    return None


# In[29]:


print_palindromes_from_list(find_DNA_palindromes(example_sequence))


# Alright. Now that that seems to work, let's try another test!

# In[30]:


new_test_file = Path("data/rosalind_revp2.txt")

(new_test_name, new_test_sequence) = read_fasta_file(new_test_file)

palindromes = find_DNA_palindromes(new_test_sequence)

print_palindromes_from_list(palindromes)

