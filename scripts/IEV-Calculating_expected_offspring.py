#!/usr/bin/env python
# coding: utf-8

# # Calculating expected offspring, from [rosalind.info](http://rosalind.info/)
# 
# (text copied from http://rosalind.info/problems/iev/)
# 
# > **Problem**  
# > For a random variable $X$ taking integer values between 1 and $n$, the expected value of $X$ is
# > $E(X) = \sum_{k=1}^{n}k×Pr(X=k)$.
# > The expected value offers us a way of taking the long-term average of a random variable over a large number of trials.
# >
# > As a motivating example, let $X$ be the number on a six-sided die.
# > Over a large number of rolls, we should expect to obtain an average of 3.5 on the die (even though it's not possible to roll a 3.5).
# > The formula for expected value confirms that
# > $E(X) = \sum_{k=1}^{6}k×Pr(X=k)=3.5$.
# >
# > More generally, a random variable for which every one of a number of equally spaced outcomes has the same probability is called a uniform random variable (in the die example, this "equal spacing" is equal to 1).
# > We can generalize our die example to find that if $X$ is a uniform random variable with minimum possible value $a$ and maximum possible value $b$, then
# > $E(X) = \frac{a+b}{2}$ 
# > You may also wish to verify that for the dice example, if $Y$ is the random variable associated with the outcome of a second die roll, then
# > $E(X+Y)=7$.
# >
# > **Given**: Six nonnegative integers, each of which does not exceed 20,000.
# > The integers correspond to the number of couples in a population possessing each genotype pairing for a given factor.
# > In order, the six given integers represent the number of couples having the following genotypes:
# >
# > 1. `AA-AA`  
# > 2. `AA-Aa`  
# > 3. `AA-aa`  
# > 4. `Aa-Aa`  
# > 5. `Aa-aa`  
# > 6. `aa-aa`
# >
# > **Return**: The expected number of offspring displaying the dominant phenotype in the next generation, 
# > under the assumption that every couple has exactly two offspring.
# >
# > **Sample Dataset**
# >
# > `1 0 0 1 0 1`
# >
# > **Sample Output**
# > `3.5`

# ## My interpretation/reasoning
# 
# 1. I am going to calculate an expected value.
# 
# 2. The calculation is a genetics exercise, involving allele frequencies.
# 
# 3. The given numbers represent pairs of different genotypes: combinations of
#     - Homozygous dominant
#     - Homozygous recessive
#     - Heterozygous
# 
# 4. The answer should be the number of offspring with the dominant phenotype.
# 
# 5. Offspring with a dominant phenotype should have at least one dominant allele.
# 
# 6. Chances of getting a dominant allele differ per parent combination:
#     1. 1 (offspring always gets a dominant allele)
#     2. 1
#     3. 1
#     4. $0.5 + 0.5 * 0.5 = 0.75 (0.5 chance to get a dominant allele from parent 1, otherwise (other 0.5) it might get a dominant allele from parent 2 with a 0.5 chance)
#     5. 0.5 (it either gets the dominant allele from parent 1 or not)
#     6. 0 (no dominant alleles, no chance of getting one)
#     
#     **Note: these chances are per offspring, and each parent pair is expected to have _two_ offspring!**
# 
# 7. Now all these chances need to be connected to their respective positions for the script to work.
# 
# Now practically, this means the code is to:
# 
#   - Open and read a text file with 6 space-separated positive integers or zero (0-20,000)
#   - These integers should be multiplied by their respective chance to produce offspring with a dominant allele (see point 6 above)
#   - Those numbers should be summed and returned as answer.

# In[1]:


def calculate_expected_offspring(numbers, offspring=2):
    """
Given a list of six positive integers and the number of offspring per pair (default=2),
calculate the expected number of offspring with a dominant allele.
    """
    dominant_offspring_likelihood = {
        1: 1 * offspring,
        2: 1 * offspring,
        3: 1 * offspring,
        4: 0.75 * offspring,
        5: 0.5 * offspring,
        6: 0 * offspring
    }
    
    #For each position in the list of 6 numbers (index 0-5),
    # multiply the number of pairs by the chance of producing
    # offspring with a dominant phenotype given the number
    # of offspring each pair produces.
    dominant_offspring = (
    numbers[0] * dominant_offspring_likelihood[1] +
    numbers[1] * dominant_offspring_likelihood[2] +
    numbers[2] * dominant_offspring_likelihood[3] +
    numbers[3] * dominant_offspring_likelihood[4] +
    numbers[4] * dominant_offspring_likelihood[5] +
    numbers[5] * dominant_offspring_likelihood[6])
    
    return(dominant_offspring)


# In[2]:


calculate_expected_offspring([1, 0, 0, 1, 0, 1])


# So far so good. Now add a function that reads numbers from a file and checks if there are six numbers between 1 and 20,000.

# In[3]:


def read_numbers(input_file):
    """
Read a file with six numbers, 
turn into a list and
confirm the length of 6.
Also check that no number exceeds 20,000.
    """
    with open(input_file, 'r') as read_file:
        for line in read_file:
            numbers_list = line.split()
            
    numbers_list = list(map(int, numbers_list))
            
    if len(numbers_list) != 6:
        print("Error, list is not length 6: %i" % len(numbers_list))
        return(None)
    else:
        pass
    
    for number in numbers_list:
        if number < 0 or number > 20000:
            print("Error! List contains a number out of range 0-20,000: %i" % number)
            return(None)
        else:
            pass
        
    return(numbers_list)


# In[4]:


read_numbers("data/Example_calculating_expected_offspring.txt")


# In[5]:


test_list = read_numbers("data/Example_calculating_expected_offspring.txt")

calculate_expected_offspring(test_list)


# Now this works too. Let's see what the function does with different numbers.

# In[6]:


calculate_expected_offspring([10, 500, 3000, 75, 15000, 938])


# As long as you provide exactly six numbers, it seems to be fine. The list check is only in the file reader, so catching errors works only when working with files. Anyway, this script seems to work.  Let's go download the dataset and solve this problem!

# In[7]:


numbers_list = read_numbers("data/rosalind_iev.txt")

calculate_expected_offspring(numbers_list)


# ## Success!!
# 
# I solved the problem. Now let's save this notebook and commit it to git.
