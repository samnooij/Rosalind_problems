#!/usr/bin/env python
# coding: utf-8

# # Enumerating k-mers lexicographically, from [Rosalind.info](https://www.rosalind.info)
# 
# (Specific exercise at: http://rosalind.info/problems/lexf/)
# 
# This exercise sounds very much like one I've done earlier: 
# Enumerating gene orders (PERM).
# 
# It seems like the only  differences are:
# 
#   - the other exercise used numbers, this one uses _letters_
#   - the other exercise required the number of possible combinations with the answer, this exercise wants _only the combinations_
#     
# So now I should be able to do basically the same thing as in [PERM](PERM-Enumerating_gene_orders.ipynb).

# In[1]:


import itertools #this library is going to do all the heavy lifting


# In[2]:


def read_alphabet_and_length(input_file):
    """
Given a path to a text file with space-separated letters on the first line,
and a number on the second line,
return the 'alphabet' as list and the number as 'length'.
    """
    first_line = True #use a trick to separate the first and second line
    
    with open(input_file, "r") as read_file:
        for line in read_file:
            if first_line:
                alphabet = line.split()
                first_line = False
            else:
                length = int(line.strip())
                
    return alphabet, length


# In[3]:


example_file = "data/Example_enumerating_k-mers_lexicographically.txt"

(example_alphabet, example_length) = read_alphabet_and_length(example_file)

print("Example alphabet: %s\nExample length: %i" %(
example_alphabet, example_length)
     )


# Alright, that works. Now let's see if I can uses the `itertools.permutations()` function to make some combinations.

# In[4]:


print(list(itertools.permutations(example_alphabet, example_length)))


# Right, that seems fine. Now I need to convert this list of tuples into separate lines of strings...

# In[5]:


for letter_combination in list(itertools.permutations(example_alphabet, example_length)):
    print("".join(list(letter_combination)))


# And with that I'm practically done, right?
# 
# Hold on..., this exercise is again not entirely clear about how the output should be ordered...
# 
# First it says the provided 'alphabet' will be ordered (which implies it may be any hypothetical order).  
# Then it says the output should be ordered alphabetically?
# 
# I will for now just go with whatever is provided and see how that works out...

# In[6]:


test_file = "data/rosalind_lexf.txt"

(test_alphabet, test_length) = read_alphabet_and_length(test_file)

for string in list(itertools.permutations(test_alphabet, test_length)):
    print("".join(list(string)))


# Hm... So this was wrong.
# 
# Let's check what the input was:

# In[7]:


print("Test alphabet: %s\nTest length: %i" %(
test_alphabet, test_length)
     )


# That looks right... What may have gone wrong then?
# 
# Oh yes, I see now that my example went wrong too. The answer I got is too short.  
# I miss all the repeats! Now how can I include those, too...?

# In[8]:


for string in list(itertools.combinations_with_replacement(test_alphabet, test_length)):
    print("".join(list(string)))


# So I used the wrong function, this one looks more like it!
# 
# Let's try the exercise again.

# In[9]:


second_test_file = "data/rosalind_lexf2.txt"

(second_test_alphabet, second_test_length) = read_alphabet_and_length(second_test_file)

for string in list(itertools.combinations_with_replacement(
    second_test_alphabet, second_test_length)):
    print("".join(list(string)))


# Oh. This is still wrong.
# 
# I thought this would be easy, but it turns out a bit more complicated...

# In[10]:


print("Second test alphabet: %s\nSecond test length: %i" %(
second_test_alphabet, second_test_length)
     )


# Oh, I see... When the 'A' has been used, it is no longer used in combinations with 'B', and so on.
# 
# So I need even more repeats than this function does by default.
# 
# (Check [this page](https://docs.python.org/3/library/itertools.html#itertools.combinations_with_replacement) for more info on the functions I tried so far.)
# 
# I will get back to this exercise and try it again later.

# Update: it looks like I need the `product()` function from itertools! (Also see: https://docs.python.org/2/library/itertools.html#itertools.product)

# In[18]:


for string in list(itertools.product(
    second_test_alphabet, repeat=second_test_length)):
    print("".join(list(string)))


# So now that I think I know how to do it, let's try it the third time:

# In[19]:


third_test_file = "data/rosalind_lexf3.txt"

(third_test_alphabet, third_test_length) = read_alphabet_and_length(third_test_file)

for string in list(itertools.product(
    third_test_alphabet, repeat=third_test_length)):
    print("".join(list(string)))

