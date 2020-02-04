#!/usr/bin/env python
# coding: utf-8

# # Mortal Fibonacci rabbits, from [rosalind.info](http://rosalind.info/)
# 
# (text copied from http://rosalind.info/problems/fibd/)
# 
# > **Problem**  
# > Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence Relations”, which followed the recurrence relation Fn=Fn−1+Fn−2
# > and assumed that each pair of rabbits reaches maturity in one month and produces a single pair of offspring (one male, one female) each subsequent month.
# >
# > Our aim is to somehow modify this recurrence relation to achieve a dynamic programming solution in the case that all rabbits die out after a fixed number of months. See Figure 4 for a depiction of a rabbit tree in which rabbits live for three months (meaning that they reproduce only twice before dying).
# >
# >   **Given**: Positive integers n≤100 and m≤20
# >
# >   **Return**: The total number of pairs of rabbits that will remain after the n-th month if all rabbits live for m months.
# >
# > Sample Dataset
# >
# > `6 3`
# >
# > Sample Output
# >
# > `4`
# 
# [Figure 4: visualisation of the example numbers](http://rosalind.info/media/problems/fibd/mortal_rabbit_tree.png)

# ## My interpretation/reasoning
# 
# 1. The _input_ will be a number `n` (number of months the 'experiment' lasts), and a number `m` (number of months after which a rabbit dies).
# 
# 2. You start with one pair of rabbits.
# 
# 3. The rabbits mature after a month, produce offspring, live another month, producing offspring again, and then they die.
#   - Each rabbit pair can produce two pairs of offspring. 
#   - Rabbits in this model have three (3) life stages: juvenile, adult 1 and adult 2. Both adult stages produce offspring. Adult 2 dies after doing so.
#   - Note that the number of life stages is equal to `m`!
#   
# Let's try a bit of code:

# In[26]:


testing = True #Am I just testing?

# Step 1: Load the input data

if testing:
    #If testing, use my own example file
    input_file = "data/Example_mortal_Fibonacci_rabbits.txt"
else:
    #Or else, use the real exercise file
    input_file = "data/Mortal_Fibonacci_rabbits.txt"

with open(input_file, 'r') as infile:
    #Read the data from the file
    input_data = infile.readline()
    
#And extract the required parameters
n = int(input_data.split()[0]) #number of months to run experiment
m = int(input_data.split()[1]) #maximum age of rabbits

print("From file %s I have read the following data:" % input_file)
print("n = %s\tm = %s\n" % (n, m))

# Step 2: Specify the model

life_stages = range(1,m+1) #for easier reading: use 1-based counting

print("Therefore, rabbits have life stages: %s" % list(life_stages))
print("and the experiment lasts for %i months." % n)

start_population = 1 #we start with one pair

print("The starting population is %i pair(s) of rabbits.\n---" % start_population)

rabbit_dict = {1: start_population} #in life_stage 1 there is now 1 rabbit pair

#Fill up the dictionary with 0 for each other stage
for stage in life_stages:
    if stage != 1:
        rabbit_dict[stage] = 0
    else:
        pass

aging_dict = {} #create a copy to move aging rabbits outside the original dictionary

print("At the start of the experiment, the rabbit dictionary is like:\n%s\n---" % rabbit_dict)

for month in range(1, n):
    print("It is now month: %i" % month)
    #Three things happen each month:
    # 1. rabbits in life stages > 1 reproduce
    newborns = 0
    for stage in rabbit_dict.keys():
        if stage > 1:
            #past stage 1 are adult, which reproduce:
            newborns += rabbit_dict[stage]
        else:
            #else they are juveniles and don't reproduce
            pass
        
    print("%i new rabbit pairs are born" % newborns)
    
    # 2. rabbits in the last life stage (m) die
    print("%i rabbit pairs die (of age)" % rabbit_dict[m])
    rabbit_dict[m] = 0
    # 3. all remaining rabbits move to the next stage
    for stage in rabbit_dict.keys():
        if stage < m:
        #Only rabbits of stages before the latest age
            aging_dict[stage + 1] = rabbit_dict[stage]
        else:
        #Others don't exist
            pass        
    
    rabbit_dict = aging_dict.copy()
    #Only add newborns to the dictionary after all the 
    # already existing rabbits have grown to the next stage.
    rabbit_dict[1] = newborns
    print("At the end of the month, we have this age distribution:\n%s\n" % rabbit_dict)
    

print("\n---\nFinally, after %i months, we end up with %i rabbit pairs" % (n, sum(rabbit_dict.values())))


# ## Implement as a function
# 
# The experiment was a success. The loose pieces of code seem to work as intended.
# 
# Now let's put all of that code into a function that does all the work:

# In[2]:


def run_experiment(n, m, start=1, testing=False):
    """
    Run the "Mortal_Fibonacci_rabbits" algorithm,
    taking as inputs a value for:
    n (int; number of months to run for)
    m (int; maximum life span of rabbits)
    start (int; number of rabbit pairs to start with;
            optional, default = 1)
    testing (bool; whether or not you are testing:
             enables printing at each step, default = False)
    """
    #From the given m, infer all the life stages:
    life_stages = range(1,m+1)
    
    #Set the starting condition in a dictionary:
    rabbit_dict = {1: start} #in life_stage 1 there is now 1 rabbit pair
    
    #Fill up the dictionary with 0 for each other stage
    for stage in life_stages:
        if stage != 1:
            rabbit_dict[stage] = 0
        else:
            pass

    aging_dict = {} #create a copy to move aging rabbits outside the original dictionary
    #This prevents errors in changing the dictionary while reading it
    
    #Now start the calculations for the experiment:
    for month in range(1, n):
        if testing:
            print("It is now month: %i" % month)
            
        #Three things happen each month:
        # 1. rabbits in life stages > 1 reproduce
        newborns = 0
        for stage in rabbit_dict.keys():
            if stage > 1:
                #past stage 1 are adult, which reproduce:
                newborns += rabbit_dict[stage]
            else:
                #else they are juveniles and don't reproduce
                pass
        
        if testing:
            print("%i new rabbit pairs are born" % newborns)

        # 2. rabbits in the last life stage (m) die
        if testing:
            print("%i rabbit pairs die (of age)" % rabbit_dict[m])
            
        rabbit_dict[m] = 0
        
        # 3. all remaining rabbits move to the next stage
        for stage in rabbit_dict.keys():
            if stage < m:
            #Only rabbits of stages before the latest age
                aging_dict[stage + 1] = rabbit_dict[stage]
            else:
            #Others don't exist
                pass        

        rabbit_dict = aging_dict.copy()
        #Only add newborns to the dictionary after all the 
        # already existing rabbits have grown to the next stage.
        rabbit_dict[1] = newborns
        
        if testing:
            print("At the end of the month, we have this age distribution:\n%s\n" % rabbit_dict)
            print("\n---\nFinally, after %i months, we end up with %i rabbit pairs" % (n, sum(rabbit_dict.values())))
    
    return(sum(rabbit_dict.values()))


# In[3]:


#Test run the new function:
print(run_experiment(6, 3, 1, True))


# ## Perform tests with other starting numbers
# 
# The function seems to work well for the numbers in the example, but how well does it work with other values for n and m?
# (Remember: the exercise states that n≤100 and m≤20.)
# 
# Let's generate some random pairs of n and m values and test the function:

# In[6]:


import random

#Run ten tests:
for test in range(10):
    n = random.randint(1,101)
    m = random.randint(1,21)
    print("n = %i, m = %i" % (n, m))
    print(run_experiment(n = n, m = m))


# Now that seems to work fine. It looks like both n and m need to be big to produce really big results. If either the experiment runs for few months (n is small), or the rabbits die young (m is small), the final result will not be very big. With high values for both (e.g. n = 56 and m = 13) the population of rabbits will grow in the billions (217,250,817,514).
# 
# ---
# 
# ## Download the dataset and do the actual exercise
# 
# I will now download the file from Rosalind.info, read the numbers into this notebook and calculate the number of rabbit pairs.

# In[7]:


#First rewrite the data import lines to a function:
def read_values_from_file(input_file, debug = True):
    """
    Read n and m values from a text file,
    assuming the file only contains two numbers
    separated by a space. E.g.:
    6 3
    
    input_file (str; path to and name of the file)
    debug (bool; print file name and values?)
    """
    with open(input_file, 'r') as infile:
        #Read the data from the file
        input_data = infile.readline()

    #And extract the required parameters
    n = int(input_data.split()[0]) #number of months to run experiment
    m = int(input_data.split()[1]) #maximum age of rabbits
    
    if debug:
        print("From file %s I have read the following data:" % input_file)
        print("n = %s\tm = %s\n" % (n, m))
        
    return(n, m)

#Then load the data by using the function
n, m = read_values_from_file("data/Mortal_Fibonacci_rabbits.txt")

#And run the experiment
print(run_experiment(n, m))

