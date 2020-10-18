#!/usr/bin/env python
# coding: utf-8

# # Enumerating gene orders, from [Rosalind.info](https://www.rosalind.info)
# 
# (Specific exercise can be found at: http://rosalind.info/problems/perm/)
# 
# ## My personal interpretation of the problem
# 
# 1. You are given a numberm which is to be the length of the sequence ($n ≤ 7$)
# 
# 2. For this number, you take all the integers up to that number (e.g. $1, 2, 3, 4, 5, 6, 7$).
# 
# 3. Then, you calculate or count the number of permutations you can make with length $n$,
# 
# 4. and finally, you print each of those permutations (on separate lines, the order does not matter).
# 
# To do that, I want my code to:
# 
# 1. read a number n between 1 and 7,
# 2. generate a list of numbers from 1 up to n,
# 3. calculate the number of permutations*,
# 4. create the list of permutations**,
# 5. and print the results.
# 
# \* : I happen to know the formula for this: permutations = (median or 'middle number' of 1 through n) * n  
# In this case that should be equal to: permutations = sum(all numbers 1 through n). E.g. for $n = 7$, that would make either $4 * 7 = 28$, or $1 + 2 + 3 + 4 + 5 + 6 + 7 = 28$.
# 
# ** : creating this list will probably be the most challenging part of this exercise. I will try and take notes in my function(s) to explain how I do it.

# In[2]:


def read_length(input_file):
    """
Read the length (n) from a file.
Assumes the length variable is in a text file with only one number.
    """
    with open(input_file, 'r') as read_file:
        for line in read_file:
            length = line
            break 
            #Only read the first line, then stop to make sure the variable
            # length is not accidentally overwritten by a second thing.
            
    #Convert length (from string) into an integer:
    length = length.strip()
    length = map(int, length)
    
    #Now make sure it is a positive integer smaller than or equal to 7:
    if length < 1 or length > 7:
        print("Length is not between 1 and 7: %i" % length)
        return None
    else:
        return length


# In[3]:


#Now let's use the 'range' function to create a list of numbers:
for number in range(1,7+1):
    print(number)
#Remember to start from 1 (instead of 0), and end with n + 1,
# because the range is non-inclusive.


# In[4]:


#To calculate the median you can use the statistics libraryb
import statistics

statistics.median(range(1,7+1))


# In[5]:


def calculate_total_permutations(length):
    """
Calculate the total number of permutations of length 'length'.
Uses the formula:
 permutations = median(1:length) * length
    """
    return statistics.median(range(1,length+1)) * length


# # In-between thoughts...
# 
# Each permutation is built up as follows:
#  - for each and any number in the list,
#  - use that number as starting position as many times as there are other numbers
#  - then for the second position, use the next number as many times as there are numbers left.
#  - continue until no numbers are left.
#  
# Is that right?
# 
# So for instance with n = 3, you get:
# 
# 1. start with `1`, which leaves 2 and 3 (2 numbers)
# 2. that means 2 lines start with `1`
# 3. the next number in line is `2`, which leaves only 3 (1 number)
# 4. so there is one `1` followed by `2`
# 5. and that leaves only 3.
# 6. Now go to the second line and pick the second number next: `3`, which leaves 2.
# 7. Add this remaining 2 to the list.
# 
# You should now have:  
# `1 2 3`  
# `1 3 2`
# 
# Repeat the steps above for each other number in the list, i.e. `2` and `3` each get to be the starting number.
# 
# Now using this process will not directly generate the numbers in the required order, so they may have to be stored in lists and/or dictionaries before printing to the terminal.

# In[6]:


def create_permutation_list(length):
    """
Create a list of each permutation of length 'length'.
    """
    numbers = range(1,length+1)
    permutations_dict = {}
    
    for number in numbers:
        permutations_dict[number] = (len(numbers) - 1) * [number]
        #For each number, reserve a list that starts with that number
        
    return permutations_dict


# In[7]:


print(create_permutation_list(3))


# In[8]:


print(list(range(1, 5)))


# In[9]:


test_list = list(range(1, 5))

print(test_list)

first_number = test_list.pop(0)

print(first_number)

print(test_list)


# ## Second thoughts
# 
# I have been searching the web for an easier solution. Turns out there are ready-made functions for this. (Of course...)  
# Also see: https://stackoverflow.com/questions/104420/how-to-generate-all-permutations-of-a-list
# 
# Probably the most useful functions are in the `itertools` library:
# ```python
# import itertools
# 
# list(itertools.permutations([1,2,3,4]))
# 
# ```
# 
# Now the only thing up to me is to print the combinations on separate lines, with spaces between the numbers.

# In[10]:


import itertools


# In[11]:


def create_permutations(length):
    """
Create the list of permutations for a given length 'length'.
This makes combinations of the numbers in the list of 1 through 'length'.
(uses 'itertools' library)
    """
    if length == 1:
        return 1
    else:
        permutations_list = itertools.permutations(range(1,length+1))
        
    return '\n'.join([' '.join(map(str, t)) for t in permutations_list])    


# In[12]:


test_tuple = (1, 2, 3)

print(test_tuple)

print(' '.join(map(str, test_tuple)))


# In[13]:


test_tuple_list = [(1, 2, 3), (2, 3, 4)]

print('\n'.join(map(str, test_tuple_list)))


# In[14]:


print(' '.join('\n'.join(map(str, test_tuple_list))))


# In[15]:


string_list = [' '.join(map(str, t)) for t in test_tuple_list]

print(string_list)

print('\n'.join(string_list))


# In[16]:


print(create_permutations(5))


# So that works, now we only need to add the number of combinations on top:

# In[20]:


len(create_permutations(3).split('\n'))


# So the final function needs to split the combinations by newlines, then count all the lines, print that number and then print all the combinations.

# In[21]:


def print_number_and_permutations(permutations):
    """
Given a newline-separated list of combinations,
return the number of combinations as well as the
original list.
    """
    number = len(permutations.split("\n"))
    return("%s\n%s" % (number, permutations))


# In[23]:


print(print_number_and_permutations(create_permutations(3)))


# Important to note: using this function requires a `print()` to get the newlines to work correctly. Otherwise "\n" will be printed as output.
# 
# Anyway, with this function I should be able to finish the assignment. Let's give it a try.

# In[25]:


test_file = "data/rosalind_perm.txt"

with open(test_file, "r") as read_file:
    for line in read_file:
        test_number = int(line.strip())
        break #read only the first line
        
number_and_permutations = print_number_and_permutations(create_permutations(test_number))

print(number_and_permutations)

