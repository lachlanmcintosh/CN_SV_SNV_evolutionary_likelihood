# import the required libraries

import numpy as np
import copy
import pickle as pkl
import re
from more_itertools import locate
import shelve
import scipy

# precomputed files
precomputed_file_folder = "/vast/scratch/users/lmcintosh/GD2/GD/"


##### STEP 1; WRITE A FUNCTION THAT CAN SIMULATE A GENOME
#####
#####
#####
#####
#####

# copied from http://www.insilicase.com/Web/Chromlen.aspx
# Percent of total (Female) genome
lengths = {0:8.18,
        1:8.04,
        2:8.60,
        3:6.33,
        4:5.98,
        5:5.65,
        6:5.25,
        7:4.84,
        8:4.64,
        9:4.48,
        10:4.45,
        11:4.38,
        12:3.78,
        13:3.52,
        14:3.32,
        15:2.94,
        16:2.61,
        17:2.52,
        18:2.11,
        19:2.07,
        20:1.55,
        21:1.64,
        22:5.13
        }

# lengths is a dictionary of rate adjustment parameters. 
# these parameters adjust the rate of how frequently SNVs occur per genome, so that longer chromosomes have a proportionally larger chance of a new SNV
# lengths is normalised so that the average rate parameter associated with each chromosome is one:
for chrom_type in lengths:
    lengths[chrom_type] *= 23/100


# count paternity takes in a list of simulated chromosomes and returns the number of each paternal type, 
# if these numbers were sorted they would be major/minor copy numbers
def count_paternity(chromosomes,paternal):
    return( len([x for x in chromosomes if paternal == x["paternal"]]) )


def simulate_single_with_poisson_timestamps_names(p_up,p_down,pre,mid,post,rate):
    # this function simulates the evolution of a cancer genome with whole chromosome copy number changes and SNV's
    # pre, mid and post model the number of epochs / cell cycles that this genome undergoes. 
    # pre is the number of epochs that occur before the first round of genome doubling, if there is one
    # mid is the number of epochs that occur after the first round of genome doubling and before the second round of genome doubling
    #   if there is no first round of genome doubling then mid is set to be -1
    # post is the number of epochs that occur after the second round of genome doubling
    #   if there is no second round of genome doubling then post is also -1
    # rate is the base rate at which SNVs occur per epoch, 
    # this rate is scaled so that each chromosome on average will gain a number of SNVs proportional to its chromosomal length

    # SNV_count will count the number of SNVs and also provide a unique number to each new SNV created
    SNV_count = 0

    # set up the initial genome
    simulated_chromosomes = {}
    for chrom_type in range(23):
        simulated_chromosomes[chrom_type] = [
                {
                    "unique_identifier":chrom_type+x,
                    "parent":-1,
                    "epoch_created":-1,
                    "paternal":x==0,
                    "SNVs":[]
                    } for x in (0,23)]
        # unique identifier is a unique number to each new chromosome created, 

    # chrom_count will count the number of chromosomes created throughout the simulation and provide each with a unique identity
    chrom_count = 46

    # now simulate forward from the initial state
    # for each epoch in the total number of epochs:
    for epoch in range(pre+mid+post+2): 
        # simulate and insert any new SNVs into the genome
        for chrom_type in simulated_chromosomes:
            for chrom in simulated_chromosomes[chrom_type]:
                # generate a random number of SNVs proportional to the length of the genome:
                additional_SNV_count = np.random.poisson(rate * lengths[chrom_type],1)[0]

                # add these SNVs to the chromosome 
                for x in range(SNV_count + 1, SNV_count + additional_SNV_count + 1):
                    chrom["SNVs"] += [{"unique_identifier":str(x),"epoch_created":epoch}] 

                # update the SNV count
                SNV_count = SNV_count + additional_SNV_count 

        # if there is a genome doubling it has to be after post, so post cannot equal -1 if mid does
        # enforce this assertion as it may be easy to forget this later on:
        assert( not(mid == -1 and post != -1) ) 

        # simulate the changes in copy number, but keep track of the SNVs
        # first see if this is a genome doubling round or an anueploidy round (they are mutually exclusive)
        if (mid != -1 and epoch == pre) or (post != -1 and epoch == pre+1+mid): 
            # if this epoch is the epoch of the first round or the second round of genome doubling
            for chrom_type in simulated_chromosomes:
                new_chromosomes = []
                for chrom in simulated_chromosomes[chrom_type]:
                    # copy each chromosome and record a unique identity and parent for it:
                    # need to deep copy or we can change the old chromosome
                    new_chromosome = copy.deepcopy(chrom) 

                    # chrom count is the next unique identifier of a chromosome
                    chrom_count += 1 
                    new_chromosome["unique_identifier"] = chrom_count
                    new_chromosome["epoch_created"] = epoch
                    new_chromosome["parent"] = chrom["unique_identifier"]
                    new_chromosomes += [new_chromosome]

                simulated_chromosomes[chrom_type] += new_chromosomes

        else:
            # this is a standard round of aneuploidy
            for chrom_type in simulated_chromosomes:
                # until a viable next epoch in evolution is simulated keep re simulating the next epoch
                # buy viable we mean that we need to have at least one copy of every chromosome
                while(True):
                    new_chromosomes = []

                    for chrom in simulated_chromosomes[chrom_type]:
                        # randomly draw from 
                        #    losing this chromosome with probability p_down, 
                        #    gaining another copy of this chromosome with probability p_up,
                        #    nothing happening with probability 1 - p_up - down:
                        change = np.random.choice([1,0,-1], 1, p=[p_up,1-p_up-p_down,p_down])

                        if change == 0:
                            # keep the old chromosome only
                            new_chromosomes += [chrom]

                        elif change == 1:
                            # create a new chromosome
                            new_chromosome = copy.deepcopy(chrom)
                            chrom_count += 1
                            new_chromosome["unique_identifier"] = chrom_count
                            new_chromosome["epoch_created"] = epoch
                            new_chromosome["parent"] = chrom["unique_identifier"]
                            new_chromosomes += [new_chromosome]
                            # but also keep the old one
                            new_chromosomes += [chrom]

                        elif change == -1: 
                            # then lose the old chromosome
                            continue

                    # ensure that there is always at least one copy of every chromosome
                    if len(new_chromosomes) != 0:
                        break 
               
                # add those chromosomes into the genome
                simulated_chromosomes[chrom_type] = new_chromosomes

    # some quick final sanity checks:
    if post == 0 or (mid == 0 and post == -1):
        # if a genome doubling round just occurred then the copy number of every chromosome will be even:
        for chrom_type in simulated_chromosomes: 
            assert(count_paternity(simulated_chromosomes[chrom_type],paternal=True) % 2 == 0)
            assert(count_paternity(simulated_chromosomes[chrom_type],paternal=False) % 2 == 0)

    for chrom_type in simulated_chromosomes:
        assert(len(simulated_chromosomes) != 0)

    return(simulated_chromosomes)

##### STEP 2; calculate the log likelihoods over the precomputed domain for total parental specific copy number
#####
#####
#####
#####
#####

# count the number of each parental specific copy number found in the genome
def count_CN_multiplicities(simulated_chromosomes):
    multiplicities = {}

    for chrom_type in simulated_chromosomes:
        CNs = [len([x for x in simulated_chromosomes[chrom_type] if paternal == x["paternal"]]) for paternal in [True,False]]

        for CN in CNs:
            if CN not in multiplicities: 
                multiplicities[CN] = 1

            else:
                multiplicities[CN] += 1

    return(multiplicities)

# for every copy number sum the precomputed values weighted against their multiplicity
# then adjust for CN’s not being able to go to zero
def CN_multiplicities_to_likelihoods(CN_multiplicities):
    # these relative references will need to be modified before publishing 
    def CN_filename(CN):
        return precomputed_file_folder + "/collated_128_p128_v3_"+str(CN)+".npy"

    def ll(copy,multiplicity):
        return np.load(CN_filename(copy)) * multiplicity

    lls = None
    for copy in CN_multiplicities:
        if lls is None:
            lls = ll(copy,CN_multiplicities[copy])

        else:
            lls += ll(copy,CN_multiplicities[copy])

    # account for the inability to lose all copies of a particular chromosome:
    likelihoods = np.exp(lls - np.log(1-np.exp(np.load(CN_filename(0)))) * CN_multiplicities[0])

    # now we want to merge these computed likelihoods with the metadata columns:
    named_likelihoods = pkl.load(open(precomputed_file_folder+"/collated_128_p128_v3_list.pickle",'rb'))

    named_likelihoods.insert(3,"likelihood",likelihoods,True)
    named_likelihoods.columns = ["p_up","p_down","path","likelihood"]
    # named_likelihoods.replace([np.inf, -np.inf], np.nan, inplace=True)
    named_likelihoods.dropna(axis=0)
    named_likelihoods.sort_values(by=['likelihood'], inplace=True, ascending=False)
    total = np.nansum(named_likelihoods["likelihood"])
    print(total)
    named_likelihoods["likelihood"] /= total
    print("best likelihoods")
    print(named_likelihoods[:][0:300].to_string())

    return(named_likelihoods)

###### STEP 3; recalculate the top branching process loglikelihoods by incorporating the SNV data under a poisson model
###### STEP 3a; calculate the SNV multiplicities of each chromosome
######
######
######
######

# first count the copy number of each SNV:
def simulated_chromosomes_to_SNV_counts(simulated_chromosomes):
    SNV_copy_counter = {}
    for chrom_type in simulated_chromosomes:
        if chrom_type not in SNV_copy_counter:
            SNV_copy_counter[chrom_type] = {}

        for chrom in simulated_chromosomes[chrom_type]:
            for SNV in chrom["SNVs"]:
                UI = SNV["unique_identifier"]

                if UI not in SNV_copy_counter[chrom_type]:
                    SNV_copy_counter[chrom_type][UI] = 1

                else:
                    SNV_copy_counter[chrom_type][UI] += 1

    return(SNV_copy_counter)


# take the copy number of each SNV and then count the number of SNVs at a particular copy number
# SNVs of a particular copy number have a good chance of having accumulated together on the same chromosome
# we will exploit this by finding the possible implied tree structures in the evolutionary history of the chromosomes
# and use a constant rate modelling assumption to allow for a best tree structure to explain the data

def SNV_counts_to_SNV_multiplicities(SNV_copy_counter):
    multiplicities = {}
    for chrom_number in SNV_copy_counter:
        if chrom_number not in multiplicities:
            multiplicities[chrom_number] = {}

        for SNV in SNV_copy_counter[chrom_number]:
            CN = SNV_copy_counter[chrom_number][SNV]

            if CN not in multiplicities[chrom_number]:
                multiplicities[chrom_number][CN] = 1

            else:
                multiplicities[chrom_number][CN] += 1

    return(multiplicities)

def count_SNV_multiplicities(simulated_chromosomes):
    return(SNV_counts_to_SNV_multiplicities(simulated_chromosomes_to_SNV_counts(simulated_chromosomes)))


##### STEP 3b; now we have SNV counts, make all possible trees that could explain those SNV counts for the given epoch structure “(pre,mid, post)”
##### 
##### 
##### 
##### 
##### 
##### To save redundancy and speed up computational time-complexity of discovering and sorting through all evolutionary 
##### trees that explain the SNV data we do not insert copy number 1’s when finding all tree structures too iterate over.
##### They can be considered implicit and inserted everywhere at the very end. 
##### Non unitary copy numbers are not necessarily 
##### The reason for this is that there is a bijection between tree structures that include copy number 1’s and 
##### have paired timeline arrays where it every copy number must have existed for non zero time except for copy number 1
##### and the tree structures that do not include copy number 1 as every leaf of every tree and force all nodes to have 
##### a non zero evolutionary history where SNVs were allowed to accumulate. The latter is easier to computationally discover.
##### The tree structures can be simply constructed in a tuple of tuples format like tree = (value, left subtree, right subtree). 

# given a collection of binary trees and a value of ‘CN’ to be inserted into these binary trees, 
# insert a node of value CN once into every possible location for each tree and return that list of all the created trees:
def insert_node(trees, CN):
    # new_trees will be a list of trees where each tree from ‘trees’ had a value of ‘CN’ inserted once somewhere within it
    new_trees = []

    for tree in trees:
        if len(tree) == 1 and CN < tree[0]: 
            # if it is a leaf node and CN is less than the value of this leaf node insert it and append it to the list of output trees:
            new_CNs = (node,tree[0]-CN)

            # make the new trees created left side heavy to allow for matching later on to the simulated truth tree, if it exists.
            new_tree = (tree[0],(max(new_CNs),),(min(new_CNs),))
            new_trees.append(new_tree)

        elif len(tree) == 3:
            # attempt to insert this value into the left and the right subtrees:
            for subtree in insert_node([tree[1]],CN): 
                new_trees.append((tree[0],subtree,tree[2]))

            for subtree in insert_node([tree[2]],CN):
                new_trees.append((tree[0],tree[1],subtree))

    return(new_trees)


# for a given tree made up of the observed SNV’s excluding copy number 1, 
# insert 1’s everywhere until 1’s are all the leaf nodes and “completed”:
def complete_tree(tree):
    if len(tree) == 1 and (tree[0] == 0 or tree[0] == 1):
        return tree

    if len(tree) == 3:
        return((tree[0],complete_tree(tree[1]),complete_tree(tree[2])))

    elif len(tree) == 2:
        return((tree[0],complete_tree(tree[1]),complete_tree([tree[0]-tree[1][0]])))

    elif len(tree) == 1:
        return((tree[0],complete_tree((tree[0]-int(tree[0]/2),)),complete_tree((int(tree[0]/2),))))
    else:
        assert(1==2) # throw a better error message than this
    return None


def complete_trees(trees):
    return([complete_tree(tree) for tree in trees])


# from the multiplicity counts of the chromosomes and the SNVs generate all the possible trees of every copy number:
def generate_trees(chrom_CNs,SNV_CNs):
    SNV_CNs.sort(reverse = True)
    chrom_CNs.sort(reverse = True)
    # print("SNV_CNs")
    # print(SNV_CNs)
    # print("chrom_CNs")
    # print(chrom_CNs)

    # initially we start with the following tree for each chromosome:
    trees = [(sum(chrom_CNs),(max(chrom_CNs),),(min(chrom_CNs),))]

    for SNV_CN in SNV_CNs:
        # to save computational complexity we generate only trees from SNVs with copy numbers greater than 1
        if SNV_CN == 1:
            continue

        # insert the node at least once in every tree
        trees_with_new_node = insert_node(trees, SNV_CN)
 
        if trees_with_new_node == []:
            # then there isn’t anywhere left to insert the new node into these trees
            # this shouldn’t happen unless this SNV_CN is also in our chrom_CNs
            assert(SNV_CN in chrom_CNs)
            continue

        if SNV_CN in chrom_CNs:
            # then it has already been inserted into the tree once in the first branch split
            # we still need to attempt inserting it again it might present in the alternative tree branch
            trees = trees + trees_with_new_node
        else:
            # then we enforce that it is inserted
            trees = trees_with_new_node

        while(True):
            # continue to reinsert this node into the tree until no new trees are found
            trees_with_node_inserted_again = insert_node(trees_with_new_node, SNV_CN)
            if trees_with_node_inserted_again == []:
                break

            trees += trees_with_node_inserted_again
            trees_with_new_node = trees_with_node_inserted_again

    # now insert the “leaf nodes” into the tree - all of which are of CN 1
    trees = complete_trees(trees)

    # find a unique set of trees:
    trees = list(set(trees))

    return(trees)

##### STEP 3c; now that all trees have been created, calculate all possible timings for each tree
##### 
##### 
##### 
##### 
##### 


def label_tree(tree, count, parents, label_to_copy):
    # what does the word count mean, need to change that to somehting more appropriate

    # change the numbers to labels "in place"
    # as we modify the tree we slowly peel of the tuple container and place a list one around it instead
    tree = list(tree)

    # we have seen count-1 nodes so far in the tree, so we can uniquely label this one as count:
    unique_label = str(count)

    # we record the copy number of this unique label for convenience later:
    label_to_copy[unique_label] = tree[0]

    # we replace the value of the node of the tree with its label instead of its CN.
    tree[0] = unique_label
 
    # now recursively label each of the left and right subtrees (if they exist):
    new_parent = unique_label

    if len(tree) >= 2:
        tree[1], count, parents, label_to_copy = label_tree(tree[1], count+1, parents, label_to_copy)
        parents[tree[1][0]] = new_parent

        if len(tree) == 3:
            tree[2], count, parents, label_to_copy = label_tree(tree[2], count+1, parents, label_to_copy)
            parents[tree[2][0]] = new_parent

    return((tree, count, parents, label_to_copy))


def get_timings_per_tree(tree,epochs,copy_dict):
    # for a given tree assign labels to it and describe the parental relationships within the tree so that we can 
    # compute the timings of nodes in arrays

    labelled_tree, count, parents, label_to_copy = label_tree(copy.deepcopy(tree),0,{},{})
    # labelled_tree simply assigns a unique number to every node in the tree
    # count is the number of nodes in the tree
    # parents is a dictionary that describe the parental relationships between the unique numbers of the nodes in the tree
    # label_to_copy is a dictionary where the keys are the unique labels and the values are the copy number of the corresponding node

    # set up an empty array to begin with, each column will represent a particular node in the tree
    timings = np.array([None]*count)
    unique_tree_labels = range(count)

    for label in unique_tree_labels:

        if label == 0:
            # then it is the root node and there is no such thing as a timing for the first bifurcation
            # it is simply a natural split do to each person inheriting one of each CN
            timings = np.tile(timings, (1,1))
            timings[:,label] = 0 

        elif label == 1 or label == 2:
            timings = np.tile(timings, (epochs,1))
            timings[:,label] = list( range(1, epochs+1) )

        else:
            parent = int( parents[ str( label ) ] ) 

            # each row in timings is a particular solution to the problem of finding bifurcation times of each node in a given tree
            for row in range(len(timings)):
                parents_time = timings[row][parent]

                # not really too sure why parents_time could be none
                assert(not (parents_time is None))

                if parents_time <= epochs and label_to_copy[str(label)] == 1: 
                    # the copy number of this node is 1 and it doesn’t bifurcate so it can exist for 0 time
                    timings_temp = np.tile(timings[row], (1,1))
                    timings_temp[:,label] = epochs 

                elif parents_time < epochs:
                    # the copy number of this node is not 1, and therefore is must bifurcate and must exist for non zero time
                    timings_temp = np.tile(timings[row], (epochs-parents_time, 1))
                    timings_temp[:,label] = list(range(parents_time+1, epochs+1))

                else:
                    continue

                # save the timings to new_timings whilst we finish filling out this column
                if row == 0:
                    new_timings = timings_temp

                else:
                    new_timings = np.vstack([new_timings, timings_temp])
            
            # this column in the timings array has been filled with times that satisfy the constraint for its respective node 
            # so now we can update the timings array and move onto the next node/column
            timings = new_timings

    return((tree, labelled_tree, count, timings, parents))

##### STEP 3d; now that all timing arrays for the nodes of each tree have been created, calculate the likelihoods from them
##### 
##### 
##### 
##### 
##### 

# this function will find the indices of the item_to_find in the list_to_check:
def find_indices(list_to_check, item_to_find):
    indices = locate(list_to_check, lambda x: x == item_to_find)
    return list(indices)


# from the perspective of calculating the SNV likelihoods we care about how long it took for branches to bifurcate.
# we can obtain these branch_lengths from get_branch_lengths
def get_branch_lengths(timings):
    tree, labelled_tree, count, timing_array, parents = timings

    branch_lengths = copy.deepcopy(timing_array)

    for child in parents:
        ch = int(child)
        p = int(parents[child])
        branch_lengths[:,ch] = branch_lengths[:,ch] - branch_lengths[:,p]
        # these branch_lengths are the lengths for each node to bifurcate

    # now we need to stack the branch lengths of the same copy numbers together:
    CNs = sorted([x for x in re.split("\(|\)|,|'", str(tree)) if x.isdigit()],reverse=True)
    unique_CNs = sorted(list(set(CNs)),reverse=True)

    for CN in unique_CNs:
        indices = find_indices(CNs,CN)
        new_stacked_branch_lengths = branch_lengths[:,indices].sum(axis=1)

        if CN == CNs[0]:
            stacked_branch_lengths = new_stacked_branch_lengths

        else:
            stacked_branch_lengths = np.vstack((stacked_branch_lengths, new_stacked_branch_lengths))

    return((CNs, lengths_array, unique_CNs, branch_lengths, np.transpose(stacked_branch_lengths)))

# WHY DOES THIS NEED TO BE TRANSPOSED?

def get_poisson_loglikelihood(lengths,counts,branch_lengths,plambda):
    A = np.log(branch_lengths.astype(float) * plambda * lengths[chrom]) * counts
    B = -np.tile( [scipy.special.gammaln(x+1) for x in counts], (branch_lengths.shape[0],1))
    C = -branch_lengths * plambda
    total = np.sum(A + B + C, axis=1)
    return(total)

def get_all_poisson_loglikelihoods_per_chr(timings,plambda,BP_probs,observed_SNV_multiplicities): # these "timings" are on a per chromosome basis
    SNV_likelihoods = []
    for i in range(len(timings)):
        CNs, lengths_array, unique_CNs, branch_lengths, stacked_branch_lengths = get_branch_lengths(timings[i])

        copies = []
        for CN in unique_CNs:
            if int(CN) not in observed_SNV_multiplicities:
                copies += [0]
            else:
                copies += [observed_SNV_multiplicities[int(CN)]]

        # put together BP and poisson likelihoods
        this_SNV_likelihood = get_poisson_loglikelihood(copies, stacked_branch_lengths, plambda) 

        this_SNV_likelihood += BP_probs[i]

        SNV_likelihoods += [this_SNV_likelihood]

    return(SNV_likelihoods)

def find_best_SNV_likelihood(plambda, timings, BP_probs):
    SNV_likelihoods = {}
    best = {}

    chroms = timings.keys()

    for chrom in chroms:
        SNV_likelihoods[chrom] = get_all_poisson_loglikelihoods_per_chr(timings[chrom], plambda, BP_probs[chrom])
        best[chrom] = (-np.Inf, 0, 0) # the second entry is the tree and the third entry is the row of that timings tree

        for tree in range(len(SNV_likelihoods[chrom])):
            the_max = max(SNV_likelihoods[chrom][tree])
            the_row = np.argmax(SNV_likelihoods[chrom][tree])

            if the_max > best[chrom][0]:
                best[chrom] = (the_max,tree,the_row)

    total = 0
    for chrom in chroms:
        total += best[chrom][0]

    return(total,best) # also need to return which tree is the best and which row of that tree is the best.    


print("START")
pre = 2
mid = 2
post = -1
p_up=0.3
p_down=0.3
rate = 10

true_pre = pre
true_mid = mid
true_post = post
true_p_up = p_up
true_p_down = p_down
true_rate = rate

simulated_chromosomes = simulate_single_with_poisson_timestamps_names(
        p_up=p_up, 
        p_down=p_down, 
        pre=pre, 
        mid=mid, 
        post=post, 
        rate=rate)

print("copynumber multiplicities")
observed_CN_multiplicities = count_CN_multiplicities(simulated_chromosomes=simulated_chromosomes)
print(observed_CN_multiplicities)

print("loglikelihoods")
likelihoods = CN_multiplicities_to_likelihoods(observed_CN_multiplicities=observed_CN_multiplicities)
print(likelihoods)

do_steps_123 = True 

def likelihoods_to_marginal_likelihoods(likelihoods):
    aggregated_likelihoods = likelihoods#.copy()
    aggregated_likelihoods["mean_p_up"]   = aggregated_likelihoods["p_up"]  *aggregated_likelihoods["likelihood"]
    aggregated_likelihoods["mean_p_down"] = aggregated_likelihoods["p_down"]*aggregated_likelihoods["likelihood"]
    aggregated_likelihoods = aggregated_likelihoods.groupby(["path"], as_index=False)[['likelihood','mean_p_up','mean_p_down']].sum()
    aggregated_likelihoods.dropna(axis=0)
    aggregated_likelihoods.sort_values(by=['likelihood'], inplace=True, ascending=False)
    total = sum(aggregated_likelihoods["likelihood"]) 
    aggregated_likelihoods["mean_p_up"] /= total
    aggregated_likelihoods["mean_p_down"] /= total
    print("best marginal likelihoods")
    print(aggregated_likelihoods[:][0:20].to_string())
    print(sum(aggregated_likelihoods["likelihood"]))

    return(aggregated_likelihood)

print("marginals")
marginals = likelihoods_to_marginal_likelihoods(likelihoods)

if do_steps_123:

    # now we need to know the indicies of these likelihoods




    # I don't think the mean values are the best values to use, but they are the easiest to compute for now
    # it would also be easy to switch to the maximum values, so perhaps do that

    # need to save the important datastructures up to hear and then just work onwards from here to speed up development
    d = shelve.open('file.txt')           # in this file you will save your variables
    #d['named_likelihoods'] = named_likelihoods             # thats all, but note the name for later.
    d['aggregated_likelihoods'] = aggregated_likelihoods
    d['simulated_chromosomes'] = simulated_chromosomes
    d['observed_CN_multiplicities'] = observed_CN_multiplicities
    d['likelihoods'] = likelihoods
    d.close()

def path_code_to_pre_mid_post(path):
    bits = [int(x) for x in path.split("G")] + [-1,-1]
    pre, mid, post = bits[0:3]
    return((pre,mid,post))

import shelve
d = shelve.open('file.txt')
#named_likelihoods = d['named_likelihoods']
aggregated_likelihoods = d['aggregated_likelihoods']
simulated_chromosomes = d['simulated_chromosomes']
observed_CN_multiplicities = d['observed_CN_multiplicities']
likelihoods = d['likelihoods']
d.close()

# STEP 4; 

print("SNV multiplicities")
observed_SNV_multiplicities = count_SNV_multiplicities(simulated_chromosomes)
print(observed_SNV_multiplicities)


path = aggregated_likelihoods["path"].iloc[0]
p_up = aggregated_likelihoods['mean_p_up'].iloc[0].round(decimals = 2)
p_down = aggregated_likelihoods['mean_p_down'].iloc[0].round(decimals = 2)
pre, mid, post = path_code_to_pre_mid_post(path)


print("calc likelihood")
print("pre: "+ str(pre))
print("mid: "+str(mid))
print("post: "+str(post))

all_trees = {}
filled_trees = {}
timings = {}
SNV_likelihoods = {}

#def count_CN_pairs(simulated_chromosomes):
#    # count how many of each chromosomal pair there are:
#    return( [len([x for x in simulated_chromosomes if x["paternal"] == paternal]) for paternal in [True,False]] )
#
#def count_CN_pairs(simulated_chromosomes):
#    chrom_CNs = {}
#    for chrom in simulated_chromosomes:
#        chrom_CNs[chrom] = count_CN_pairs(simulated_chromosomes=simulated_chromosomes)
#
#    return(chrom_CNs)

chrom_CNs = count_chrom_CN_multiplicities(simulated_chromosomes)
for chrom in observed_SNV_multiplicities:
    print("###CHROM: "+str(chrom))

    print("chromosomal CNs")
    print(chrom_CNs[chrom])
    print("observed_SNV_multiplicities")
    print(observed_SNV_multiplicities[chrom])
    all_trees[chrom] = generate_trees(
            chrom_CNs= list(chrom_CNs[chrom].values()),
            SNV_CNs= list(observed_SNV_multiplicities[chrom].keys())
            )
    print("trees:")
    for tree in all_trees[chrom]:
        print(tree)

    epochs = pre+mid+post+(mid>=0)+(post>=0)

    timings[chrom] = [get_timings_per_tree(x,epochs,observed_SNV_multiplicities[chrom]) for x in all_trees[chrom]]

    print("timings")
    print(timings[chrom])

    print("num trees: "+ str(len(all_trees[chrom])))
    print("num filled trees: " + str(len(all_trees[chrom])))


print(timings.keys())
def objective_function_SNV_loglik(plambda,timings,BP_probs):
    total,best = find_best_SNV_likelihood(plambda,timings,BP_probs)
    print(best)
    return(-total)


# from scipy.optimize import NonlinearConstraint, Bounds
# output = scipy.optimize.minimize(fun=objective_function_SNV_loglik,x0=5,args=(timings))
# output = scipy.optimize.differential_evolution(func=objective_function_SNV_loglik,bounds=Bounds([0.], [20.]),args=(timings))
# need to reintroduce some sort of non array search for both the rate parameter AND p_up AND p_down


print(timings)

def get_path_code(code_list):
    output = ""
    count = 0
    for i in range(len(code_list)):
        if code_list[i] == "A":
            count += 1
        if code_list[i] == "GD":
            output += str(count)
            count = 0
            output += "G"
    output += str(count)
    return(output)

def timing_struct_to_BP_likelihood(data, timings, chrom, pre, mid, post, up, down):

    all_BP_probs = []

    for these_timings in timings[chrom]:
        labels, lengths_array, unique_labels, branch_lengths, stacked_branch_lengths = get_branch_lengths(these_timings)

        path = []
        if pre > -1:
            path += ["A"]*pre
        if mid > -1:
            path += ["GD"]
        if mid > 0:
            path += ["A"]*mid
        if post > -1:
            path += ["GD"]
        if post > 0:
            path += ["A"]*post

        print("timings")
        print(these_timings)
        print("labels")
        print(labels)
        print("lengths")
        print(lengths_array)
        #>>> data = pickle.load(open("pre_mat129_u65_d10.precomputed.pickle",'rb'))
        # in the keys are all the possible paths...
        # each of these is a matrix that you can use to calculate the possible paths
        
        ends = these_timings[3]
        starts = ends - lengths_array

        paths = np.zeros(ends.shape, dtype=float, order='C')
        for row in range(lengths_array.shape[0]):
            for col in range(lengths_array.shape[1]):
                #print(path)
                #print((row,col))
                these_paths = path[starts[row][col]:ends[row][col]]
                #print(these_paths)
                path_code = get_path_code(these_paths)
                #print(path_code)

                if labels[col] == '1':
                    likelihood = data[path_code][1][1]

                else:
                    likelihood = data[path_code][1][2]

                paths[row][col] = likelihood

        print("starts")
        print(starts)
        print("ends")
        print(ends)
        print("codes")
        print(paths)

        ll = np.log(paths)
        print(ll)
        BP_probs = np.sum(ll, axis=1)
        print(BP_probs)
        all_BP_probs += [BP_probs]

    return(all_BP_probs)


data = pkl.load(open("../precomputed/store_pre_pickle/pre_mat129_u"+str(int(100*p_up))+"_d"+str(int(100*p_down))+".precomputed.pickle",'rb'))
print("###########\n"*10)
BP_probs = {}
for chrom in timings.keys():
    print(chrom)
    BP_probs[chrom]  = timing_struct_to_BP_likelihood(data,timings,chrom,pre=pre,mid=mid,post=post,up=p_up,down=p_down)



# it really looks like the likelihood contributed from the BP is much lower than the ones contributed by from the SNVs
outputs = []
for plambda in range(100):
    lam = plambda/5
    outputs += [(lam,objective_function_SNV_loglik(lam,timings,BP_probs))]


#iterating though linearly doesn't seem to be too bad. It isn't fast but it isn't too bad.
# maybe we can add in the BP likelihoods now. 
print(outputs.sort(key=lambda x: x[1]))
print(outputs)


# now we need to modify it so that we are also printing out the best subtree for each tree
result = find_best_SNV_likelihood(outputs[0][0],timings,BP_probs)
print(result)

def tidy_tree_timings(tree,parents):
    return(tree)


for chrom in result[1].keys():
    tree,labelled_tree,count,timing_array,parents = timings[chrom][0]
    print("##### chrom: " + str(chrom))
    print("loglik: " + str(result[1][chrom][0]))
    these_timings = timings[chrom][result[1][chrom][1]]
    print("tree: " + str(these_timings[0]))
    print("tree timings est: " + str(tidy_tree_timings(these_timings[3][result[1][chrom][2]],parents)))
    print("tree timings truth: " + str(simulated_chromosomes[chrom]))
    print("SNV timings: " + str(these_timings[-1]))
    # i guess we just want to get the t_chr stuff out...
    # simulation[2] tells you when each chromosome was created... but we need to know its type
    # simulation[0] has the ssame len as simulation[2]
    # this is great, but now we need to record whom is the parent of whom in the simulation. can't remember how to pull it out from the snv list. that might help but better to bake it in from the start
    # then I want to compare the trees and the timings. 
