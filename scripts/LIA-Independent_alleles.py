#!/usr/bin/env python
# coding: utf-8

# # Independent alleles, from [rosalind.info](http://rosalind.info/)
# 
# (text copied from http://rosalind.info/problems/lia/)
# 
# > **Mendel's Second Law**
# >
# > Recall that Mendel's first law states that for any factor, an individual randomly assigns one of its two alleles to its offspring.
# > Yet this law does not state anything regarding the relationship with which alleles for different factors will be inherited.
# >
# > After recording the results of crossing thousands of pea plants for seven years, Mendel surmised that alleles for different factors are inherited with no dependence on each other.
# > This statement has become his second law, also known as the law of independent assortment.
# >
# > What does it mean for factors to be "assorted independently?"
# > If we cross two organisms, then a shortened form of independent assortment states that if we look only at organisms having the same alleles for one factor, then the inheritance of another factor should not change.
# >
# > For example, Mendel's first law states that if we cross two $Aa$ organisms, then 1/4 of their offspring will be $aa$, 1/4 will be $AA$, and 1/2 will be $Aa$.
# > Now, say that we cross plants that are both heterozygous for two factors, so that both of their genotypes may be written as $Aa Bb$.
# > Next, examine only $Bb$ offspring: Mendel's second law states that the same proportions of $AA$, $Aa$, and $aa$ individuals will be observed in these offspring.
# > The same fact holds for $BB$ and $bb$ offspring.
# >
# > As a result, independence will allow us to say that the probability of an $aa BB$ offspring is simply equal to the probability of an $aa$ offspring times the probability of a $BB$ organism, i.e., 1/16.
# >
# > Because of independence, we can also extend the idea of Punnett squares to multiple factors, as shown in Figure 1. We now wish to quantify Mendel's notion of independence using probability.
# >
# > [Figure 1](http://rosalind.info/media/problems/lia/dihybrid_cross.png).
# > Mendel's second law dictates that every one of the 16 possible assignments of parental alleles is equally likely.
# > The Punnett square for two factors therefore places each of these assignments in a cell of a 4 X 4 table.
# > The probability of an offspring's genome is equal to the number of times it appears in the table, divided by 16.
# 
# > **Problem**
# >
# > Two events $A$ and $B$ are independent if $Pr(A and B)$ is equal to $Pr(A)×Pr(B)$.
# > In other words, the events do not influence each other, so that we may simply calculate each of the individual probabilities separately and then multiply.
# >
# > More generally, random variables $X$ and $Y$ are independent if whenever $A$ and $B$ are respective events for $X$ and $Y$, $A$ and $B$ are independent (i.e., $Pr(A and B)=Pr(A)×Pr(B)$).
# >
# > As an example of how helpful independence can be for calculating probabilities, let $X$ and $Y$ represent the numbers showing on two six-sided dice.
# > Intuitively, the number of pips showing on one die should not affect the number showing on the other die.
# > If we want to find the probability that $X+Y$ is odd, then we don't need to draw a tree diagram and consider all possibilities.
# > We simply first note that for $X+Y$ to be odd, either $X$ is even and $Y$ is odd or $X$ is odd and $Y$ is even.
# > In terms of probability, $Pr(X+Y is odd) = Pr(X is even and Y is odd) + Pr(X is odd and Y is even)$.
# > Using independence, this becomes $[Pr(X is even) × Pr(Y is odd)] + [Pr(X is odd) × Pr(Y is even)]$, or $(\frac{1}{2})^2 + (\frac{1}{2})^2 = \frac{1}{2}$.
# > You can verify this result in Figure 2, which shows all 36 outcomes for rolling two dice.
# >
# > [Figure 2](http://rosalind.info/media/problems/lia/two_dice.png). The probability of each outcome for the sum of the values on two rolled dice (black and white), broken down depending on the number of pips showing on each die.
# > You can verify that 18 of the 36 equally probable possibilities result in an odd sum.
# >
# > **Given**: Two positive integers $k (k≤7)$ and $N (N≤2k)$.
# > In this problem, we begin with Tom, who in the 0th generation has genotype Aa Bb.
# > Tom has two children in the 1st generation, each of whom has two children, and so on.
# > Each organism always mates with an organism having genotype Aa Bb.
# >
# > **Return**: The probability that at least $N$ Aa Bb organisms will belong to the $k$-th generation of Tom's family tree (don't count the Aa Bb mates at each level). 
# > Assume that Mendel's second law holds for the factors.
# >
# > **Sample Dataset**
# > ```
# > 2 1
# > ```
# > **Sample Output**
# > ```
# > 0.684
# > ```

# ## My interpretation/reasoning
# 
# 1. This is an exercise of genetics and probabilities.
# 
# 2. The exercise focuses on double heterozygous organisms ($AaBb$).
# 
# 3. We start from one individual, Tom, with genotype AaBb.
# 
# 4. Tom mates with an AaBb and has two children in generation 1 (k = 1).
# 
# 5. In the next generation, those children each mate with an AaBb and have two children of their own, until generation k (<=7).
# 
# 6. For a given number N, I am going to calculate the probability that at least that many organisms have genotype AaBb. (Assuming independent inheritence.)
# 
# 7. The example reads: After two (k = 2) generations, there will be 2 * 2 = 4 organisms. The probability that at least 1 (N = 1) has genotype AaBb is 0.684.
# 
# 8. For the first generation, there are two organisms. Their genotypes can be any of the 16 shown in [figure 1](http://rosalind.info/media/problems/lia/dihybrid_cross.png). 4/16 possibilities (1/4) are AaBb. Then for next generations, probabilities have to be calculated for each possible genotype that mate with an AaBb. It is probably best to calculate all and finally return only the probability for AaBb in the script.
# 
# 9. Since we assume independent inheritence, we may calculate probabilities for Aa and Bb separately and multiply in the end.
# 
# 10. Chances for Aa and Bb should be identical given the situation, so only one needs to be calculated and can be squared in the end.
# 
# 11. The difficulty for me now seems to in the 'at least N organisms have genotype AaBb' part. 
#     - The number of organisms doubles per generation.
#     - Probabilities are equal for all organisms
#     - The probability for 'at least N organisms' should be the sum of probabilities for N organisms in that generation.
#     - N cannot be larger than 2 * k, because that is the total number of organisms.
#     
# The difficulty seems to be mostly in the fact that organisms in generation 2 and up may have a parent from 'the Tom line' that was not AaBb. In generation 2, Aa may be the offspring of AA x Aa, Aa x Aa or aa x Aa. All of these are viable possibilities. Then for generation 3, the 'path of ancestors' becomes even more complicated, because in the generation before organisms may have had any genotype.
# 
# However, this may not really matter in the end. For each generation, there are three possible genotypes: AA, Aa and aa. Each with their respective probability.
# In the first generation, probabilities are: AA = 1/4, Aa = 1/2, aa = 1/4. Then for the second generation, the probabilities are like this:
# 
# AA = 1/2 * (AA first) + 1/4 * (Aa first) + 0 * (aa first). Where the 'first' refers to the probabilities observed in the previous generation. This same principle should apply to each generation and each genotype. Because the 'mate' has a set genotype (heterozygous), the chances for each genotype are always the same. From the second generation onwards, the only added difficulty is summing the probabilities from all 3 previous genotypes. Let's see how that can work in code...
# 
# Now practically, this means the code is to:
# 
#   - Open and read a text file with two numbers: k (generations) and N (number of individuals with genotype AaBb)
#   - For the first generation, the probabilities per organism are known: AA = 1/4, Aa = 1/2, aa = 1/4
#   - For each next generation, the probabilities are calculated based on previous probabilities (so the probabilities have to be stored, at least temporarily)
#   - This calculation has to be repeated k times.
#   - The final answer will be the probability of Aa (which is stored) at timepoint k squared (to incorporate Bb for which probabilities should be identical), to the power N (for at least N organisms have genotype Aa) PLUS that same chance for each value greater than N up to 2 * k!! (**at least N means including any number greater than N!**)

# In[101]:


def read_numbers(input_file):
    """
Read numbers k (generations) and N (number of organisms) from a file.
Confirm that k is no greater than 7 and N is no greater than 2 * k.
    """
    with open(input_file, 'r') as read_file:
        for line in read_file:
            numbers = line.split()
            generation = int(numbers[0])
            number = int(numbers[1])
            
    if generation > 7 or generation < 1:
        print("Generation is outside range 1-7: %i" % generation)
        return(None)
    else:
        pass
    
    if number > (generation ** 2) or number < 0:
        print("Number 'N' is too great, or below zero: %i" % number)
        return(None)
    else:
        pass
    
    return(generation, number)

def calculate_probabilities(previous_probabilities):
    """
Calculate the probabilities for each genotype (AA, Aa, aa)
based on the probabilities of the previous generation.
    """
    probability_AA = (previous_probabilities["AA"] * 1/2 +
                      previous_probabilities["Aa"] * 1/4 +
                      previous_probabilities["aa"] * 0)
    
    probability_Aa = (previous_probabilities["AA"] * 1/2 +
                      previous_probabilities["Aa"] * 1/2 +
                      previous_probabilities["aa"] * 1/2)
            
    probability_aa = (previous_probabilities["AA"] * 0 +
                      previous_probabilities["Aa"] * 1/4 +
                      previous_probabilities["AA"] * 1/2)
            
    new_probabilities = { "AA": probability_AA,
                          "Aa": probability_Aa,
                          "aa": probability_aa }
    
    return(new_probabilities)

def solve_independent_alleles(input_file):
    """
Function that reads an input file with numbers for k (generations) 
and N (number of organisms) to solve the exercise and return
the probability that at least N organisms in generation k have
the heterozygous genotype Aa.
    """
    (generation, number) = read_numbers(input_file)
    
    #First calculate the probabilities for the 3 genotypes in generation k
    initial_probabilities = { "AA": 1/4, "Aa": 1/2, "aa": 1/4 }
    
    if generation == 1:
        probabilities = initial_probabilities
        
        probability_Aa = probabilities["Aa"]
        
        print("The probability that one organism has genotype Aa is: %f" % probability_Aa)
        
    elif generation > 1:
        probabilities = initial_probabilities
        
        for n in range(2, generation):
            #From the second generation onwards,
            # calculate the probabilities for each genotype.
            new_probabilities = calculate_probabilities(probabilities)
            #Overwrite the 'probabilities' variable so that it is recalculated
            # for each step in this loop (number of generations from 2).
            probabilities = new_probabilities
            
            print(new_probabilities)
        
    else:
        return("Generation seems to be a number < 1: %i" % generation)
        
    if number < generation * 2:
        #If there are multiple solutions to 'at least N in k*2', e.g.
        # if N = 2 and k = 4, there are 8 organisms. At least 2/8 means it can
        # be 3, 4, 5 ... 8 as well, so add those probabilities.
        total_probability = 0

        for n in range(number, generation * 2):
            probability = probabilities["Aa"] ** 2 ** n
            #The probability is the probability for Aa squared,
            # to incorporate the identical probability of Bb,
            # then times the number of organisms that have to have
            # that exact genotype.
            total_probability += probability
            
            print(total_probability)

    else:
        #N is equal to 2*k, i.e. there is only one solution: N is all organisms in generation k.
        total_probability = probabilities["Aa"] ** 2 ** number
        #The probability is the probability for Aa squared,
        # to incorporate the identical probability of Bb,
        # then times the number of organisms that have to have
        # that exact genotype.

    return(total_probability)


# In[10]:


print(solve_independent_alleles("data/Example_independent_alleles.txt"))


# It looks like I'm underestimating the probabilities. Something is going wrong. Does it have to do with possible combinations of organisms? That if at least one out of two has to have that genotype, that is may be number one, number two, or both (a total of three possibilities)? I don't know. I will look at it later.

# These combinations can be calculated with the  `comb` function from the `math` module. Also see https://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python for more information on the topic. I think this is what I need to incorporate into my function to solve the problem.

# In[1]:


from math import comb


# In[21]:


def solve_independent_alleles2(input_file):
    """
    -UPDATED-
Function that reads an input file with numbers for k (generations) 
and N (number of organisms) to solve the exercise and return
the probability that at least N organisms in generation k have
the heterozygous genotype Aa.
    """
    (generation, number) = read_numbers(input_file)
    
    #First calculate the probabilities for the 3 genotypes in generation k
    initial_probabilities = { "AA": 1/4, "Aa": 1/2, "aa": 1/4 }
    
    if generation == 1:
        probabilities = initial_probabilities
        
        probability_Aa = probabilities["Aa"]
        
        print("The probability that one organism has genotype Aa is: %f" % probability_Aa)
        
    elif generation > 1:
        probabilities = initial_probabilities
        
        for n in range(1, generation):
            #From the second generation onwards,
            # calculate the probabilities for each genotype.
            new_probabilities = calculate_probabilities(probabilities)
            #Overwrite the 'probabilities' variable so that it is recalculated
            # for each step in this loop (number of generations from 2).
            probabilities = new_probabilities
            
            print(new_probabilities)
        
    else:
        return("Generation seems to be a number < 1: %i" % generation)
        
    if number < generation * 2:
        #If there are multiple solutions to 'at least N in k*2', e.g.
        # if N = 2 and k = 4, there are 8 organisms. At least 2/8 means it can
        # be 3, 4, 5 ... 8 as well, so add those probabilities.
        total_probability = 0

        for n in range(number, generation * 2):
            print(probabilities)
            probability = (probabilities["Aa"] * comb(generation * 2, n)) ** 2
            #The probability is the probability for Aa squared,
            # to incorporate the identical probability of Bb,
            # then times the number of organisms that have to have
            # that exact genotype.
            total_probability += probability
            
            print(total_probability)

    else:
        #N is equal to 2*k, i.e. there is only one solution: N is all organisms in generation k.
        total_probability = probabilities["Aa"] ** 2 ** number * comb(generation * 2, number)
        #The probability is the probability for Aa squared,
        # to incorporate the identical probability of Bb,
        # then times the number of organisms that have to have
        # that exact genotype.

    return(total_probability)


# In[22]:


print(solve_independent_alleles2("data/Example_independent_alleles.txt"))


# In[8]:


1.390625 / 2


# It seems that the chances of each genotype are identical each generation. The most important thing to take into account then is that when at least 1 organism in generation 2 is to be Aa(Bb), then that can happen in 4 possible ways (like 1000, 0100, 0010, 0001). AND it also counts if 2, 3 or 4 organisms have that genotype, so the chances for each of those outcomes should be added.

# In[25]:


chance_of_1 = 0.5 ** 2
combinations_of_1 = comb(4, 1)

total_probability_1 = chance_of_1 * combinations_of_1

print(total_probability_1)


# But this is too simple. This does not take into account that when the 1 does have genotype Aa, the other three do NOT. (Right?)

# In[26]:


chance_of_1 = 0.5 ** 2 * (3 * 0.5 ** 2)
combinations_of_1 = comb(4, 1)

total_probability_1 = chance_of_1 * combinations_of_1

print(total_probability_1)


# In[29]:


1/4 * 3/4 ** 3 * comb(4,1) + 1/4 ** 2 * 3/4 ** 2 * comb(4,2) + 1/4 ** 3 * 3/4 ** 1 * comb(4,3) + 1/4 ** 4


# In[30]:


total = 4
number = 1

chance_of_1 = 1/4
chance_not = 3/4
combinations_of_1 = comb(total, number)

print(combinations_of_1)


# In[31]:


1/4 * (3/4) ** 3 * 4


# In[32]:


27/64


# In[33]:


1/16 * 9/16 * 6


# In[34]:


1/4 ** 3 * 3/4 * 4


# In[50]:


1/4 ** 4 * (3/4) ** 0 * 1


# In[38]:


0.421875 + 0.2109375 + 0.046875 + 0.00390625


# In[44]:


3/4 ** 0


# In[45]:


from math import pow


# In[48]:


pow(3/4, 0)


# In[49]:


(3/4) ** 0


# With some experimentation and calculating on paper, I have come to the conclusion that:
#  - chances for each genotype remain the same over generations (so do not need to be recalculated)
#  - with numbers k and N it is actually relatively simple to calculate the chance of 'at least N organisms have a heterozygous genotype'
# 
# The calculation is as follows:
# The number of organisms per generation (k) = 2 * k
# The number of organisms that match the 'at least N' criterion is a list ranging from N  to 2 * k - for each of these outcomes you need to calculate the probability and sum all these probabilities
# For each of these numbers, there can be some combinations of organisms that give this outcome, e.g. when exactly one organism is heterozygous in a population of four, this can be represented as 1000, 0100, 0010 and 0001. So 4 possibilities: multiply the probability by 4.
# In each of those situations, the assumption is that the other organisms are NOT heterozygous.
# The chance of being heterozygous is 1/4, the chance of being not heterozygous is, therefore, 3/4. (Or 1 - 1/4)
# We then get: probability(heterozygous) to the power [number of organisms] times probability(not heterozygous) to the power [number of organisms] times the number of combinations
# For example in case of 1 organism out of a total of 4 (generation 2):
#   1/4 ** 1 * 3/4 ** 3 * 4 = 27/64 or 0.421875
# This is the chance that exactly 1 organism is heterozygous. Then add the chances for exactly 2, 3 or 4 (N+1 - 2 * k), and that should give the right answer!
# 
# Now to put this into a function...

# In[81]:


def solve_independent_alleles3(generation, number):
    """
    -UPDATED AGAIN-
Function that reads an input file with numbers for k (generations) 
and N (number of organisms) to solve the exercise and return
the probability that at least N organisms in generation k have
the heterozygous genotype Aa.
    """
    #The probabilities for each genotype are the same throughout generations,
    # because each organisms mates with a heterozygous organism.
    probabilities = { "AA": 1/4, "Aa": 1/2, "aa": 1/4 }
    
    if number < generation * 2:
        #If there are multiple solutions to 'at least N in k*2', e.g.
        # if N = 2 and k = 4, there are 8 organisms. At least 2/8 means it can
        # be 3, 4, 5 ... 8 as well, so add those probabilities.
        total_probability = 0

        for n in range(number, generation * 2 + 1):
            probability = (probabilities["Aa"] ** 2) ** n * (1 - probabilities["Aa"] ** 2) ** (2 * generation - n) * comb(generation * 2, n)
            #The probability is the probability for Aa squared,
            # to incorporate the identical probability of Bb,
            # then to the power of the number of organisms that have to have
            # that exact genotype (n). Then this is multiplied by the chance
            # that an organism has another genotype (1 - probabilities["Aa"] ** 2),
            # to the power of the number of organisms that have this genotype
            # (2 * generation - n) = (2k - n).
            # And finally, the outcome is multiplied by the number of combinations
            # that are possible for the given scenario.
            total_probability += probability

    else:
        #N is equal to 2*k, i.e. there is only one solution: N is all organisms in generation k.
        total_probability = (probabilities["Aa"] ** 2) ** number * (1 - probabilities["Aa"] ** 2) ** 0 * comb(generation * 2, number)
        #The probability is the probability for Aa squared,
        # to incorporate the identical probability of Bb,
        # then times the number of organisms that have to have
        # that exact genotype.

    return(total_probability)


# In[84]:


(generation, number) = read_numbers("data/Example_independent_alleles.txt")
solve_independent_alleles3(generation, number)


# In[90]:


solve_independent_alleles3(7, 7)


# For some reason this feels like I have created a very roundabout way of saying something like this formula:
# 
# $E(X) = \sum_{k=N}^{2g} (\frac{1}{4})^N × (\frac{3}{4}) ^ {2g - N} × \binom{N}{2g}$
# 
# where 'g' = the generation (k in the example)
# 
# Now I am not entirely sure if this mathematical notation is correct, but let's try writing this as a 'simpler' function:

# In[103]:


def solve_independent_alleles_short(k, N):
    """
The same formula as above, written shorter.
    """
    return(sum([(1/4) ** n * (3/4) ** (2*k - n) * comb(2*k, n) for n in range(N, 2*k + 1)]))


# In[96]:


solve_independent_alleles_short(2, 1)


# In[97]:


round(solve_independent_alleles_short(2, 1), 3)


# Just in case it needs to be rounded to 3 decimals, I can use the round function.
# 
# ---
# 
# Now that this seems to work, let's try the real exercise!

# In[104]:


(generation, number) = read_numbers("data/rosalind_lia.txt")

solve_independent_alleles_short(generation, number)


# This clearly didn't work... I think I got it wrong with the number of organisms per generation. I made it 2 * generation, but it should be $generation ^ 2$.
# 
# No! Perhaps I mixed up the number, put them the wrong way round? That also doesn't seem to be the case. First generation, then number.
# 
# Now what could have gone wrong...?

# I just spotted a major overlook in the number of organisms per generation! I used the calculation: organisms = 2 * k, which coincidentally works for generations 1 and 2 (which works well enough with the example).
# However, the number of organisms should be $2^k$. So the function then becomes:

# In[105]:


def solve_independent_alleles_short2(k, N):
    """
The same formula as above, written shorter.
    """
    return(sum([(1/4) ** n * (3/4) ** (2**k - n) * comb(2**k, n) for n in range(N, 2**k + 1)]))


# And the formula more like:
# 
# $E(X) = \sum_{k=N}^{2^g} (\frac{1}{4})^N × (\frac{3}{4}) ^ {2^g - N} × \binom{N}{2^g}$

# In[106]:


(generation, number) = read_numbers("data/rosalind_lia.txt")

solve_independent_alleles_short2(generation, number)


# Now I think that has solved the issue I had. Let's see if I can pass the test now:

# In[107]:


(generation, number) = read_numbers("data/rosalind_lia2.txt")

solve_independent_alleles_short2(generation, number)


# ## Success!!
# 
# Yes, that worked. Now let's save all of this and commit to git.
