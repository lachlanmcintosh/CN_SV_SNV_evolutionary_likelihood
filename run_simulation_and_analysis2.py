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

##### STEP 1b; from the simulated genome create a tree
#####
#####
#####
#####
#####

# now make a structure to compare the truth tree to the found tree
def insert_node_into_truth_tree(tree,node):
    print("node")
    print("\t"+str(node))
    print("tree")
    print("\t"+str(tree))
    assert(node["unique_identifier"] != tree["unique_identifier"])

    if node["parent"] == tree["unique_identifier"]:
        if tree["child"] == None:
            tree["complement"] = copy.deepcopy(tree)
            tree["complement"]["epoch_created"] = node["epoch_created"]
            # because tree is already a copoy of itself if should already have child and complement set to None

            tree["child"] = copy.deepcopy(node)
            tree["child"]["child"] = None
            tree["child"]["complement"] = None

        else:
            tree["complement"] = insert_node_into_truth_tree(tree["complement"],node)

    else:
        # insert it below child or complement (but not at that level)
        if tree["child"] != None:
            tree["child"] = insert_node_into_truth_tree(tree["child"],node)

        if tree["complement"] != None:
            tree["complement"] = insert_node_into_truth_tree(tree["complement"],node)

    return(tree)


# now that the truth tree is created for each chromosomes,
#   remove the SNV lists themselves, 
#   count the copynumber of each node and how many unique SNVs there are at that copy number.

# create a recursive function to insert the correct copy number at each node in the tree
def add_copynumber_to_tree(tree):
    if tree["child"] == None:
        assert(tree["complement"] == None)
        tree["copy_number"] = 1

    else:
        tree["child"] = add_copynumber_to_tree(tree["child"])
        tree["complement"] = add_copynumber_to_tree(tree["complement"])
        tree["copy_number"] = tree["child"]["copy_number"] + tree["complement"]["copy_number"]

    return(tree)

# create a recursive function to find the correct number of SNVs at a particular copy at each node in the tree:
def add_SNV_multiplicity_to_tree(tree):
    if tree["child"] == None:
        assert(tree["complement"] == None)

        count = 0
        for SNV in tree["SNVs"]:
            if SNV["epoch_created"] >= tree["epoch_created"]:
                count += 1

        tree["SNV_multiplicity"] = count

    else:
        assert(tree["child"]["epoch_created"] == tree["complement"]["epoch_created"])

        # first fill out each branch of the tree:
        tree["child"] = add_SNV_multiplicity_to_tree(tree["child"])
        tree["complement"] = add_SNV_multiplicity_to_tree(tree["complement"])
        
        # now fill out this node:
        # (we can use the "epoch_created" tag on each SNV to calculate this)
        count = 0
        for SNV in tree["SNVs"]:
            if SNV["epoch_created"] >= tree["epoch_created"] and SNV["epoch_created"] < tree["child"]["epoch_created"]:
                count += 1

        tree["SNV_multiplicity"] = count 

    return(tree)


def remove_SNVs_from_tree(tree):
    tree.pop("SNVs")

    if tree["child"] != None:
        assert(tree["complement"] != None)
        tree["child"] = remove_SNVs_from_tree(tree["child"])
        tree["complement"] = remove_SNVs_from_tree(tree["complement"])

    return(tree)

def create_truth_trees(simulated_chromosomes):
    #there is a tree fro every chromosome
    trees = {}
    for chrom_type in simulated_chromosomes:
        # first sort the nodes by the order they need to be inserted in:
        sorted_list = sorted([(x["unique_identifier"],x) for x in simulated_chromosomes[chrom_type]])

        # create the root node of the tree for this chrom_type
        tree = {'unique_identifier':-1,
                'parent':None,
                'epoch_created':None,
                'paternal':None,
                'child':None,
                'complement':None,
                'SNVs':[]
                } 

        # insert all nodes and add metadat to tree and simplify:
        for new_node in sorted_list:
            trees[chrom_type] = insert_node_into_truth_tree(tree,new_node[1])
           
        #print(trees[chrom_type])
        trees[chrom_type] = add_copynumber_to_tree(trees[chrom_type])
        #print(trees[chrom_type])
        trees[chrom_type] = add_SNV_multiplicity_to_tree(trees[chrom_type])
        #print(trees[chrom_type])
        trees[chrom_type] = remove_SNVs_from_tree(trees[chrom_type])
        #print(trees[chrom_type])

    return(trees)


def CN_tree_from_truth_tree(truth_tree):

    if truth_tree["child"] != None and truth_tree["complement"] != None:
        child_tree = CN_tree_from_truth_tree(truth_tree["child"])
        complement_tree = CN_tree_from_truth_tree(truth_tree["complement"])
        CN_tree = [truth_tree["copy_number"], child_tree,complement_tree] 

    elif truth_tree["child"] != None:
        child_tree = CN_tree_from_truth_tree(truth_tree["child"])
        CN_tree = [truth_tree["copy_number"], child_tree] 

    elif truth_tree["complement"] != None:
        complement_tree = CN_tree_from_truth_tree(truth_tree["complement"])
        CN_tree = [truth_tree["copy_number"], complement_tree] 

    else:
        CN_tree = [truth_tree["copy_number"]]

    return(CN_tree)


def make_left_heavy(tree):
    assert(len(tree) == 1 or len(tree) == 3)

    if len(tree) == 1:
        return(tree)
    else:
        if tree[1][0] < tree[2][0]:
            return([tree[0],make_left_heavy(tree[2]),make_left_heavy(tree[1])])
        else:
            return([tree[0],make_left_heavy(tree[1]),make_left_heavy(tree[2])])


def CN_trees_from_truth_trees(truth_trees):
    for chrom_type in truth_trees:
        truth_trees[chrom_type] = CN_tree_from_truth_tree(truth_trees[chrom_type])
        truth_trees[chrom_type] = make_left_heavy(truth_trees[chrom_type])
    return(truth_trees)


##### STEP 2; calculate the log likelihoods over the precomputed domain for total parental specific copy number
#####
#####
#####
#####
#####

# count the number of each parental specific copy number found in the genome
def count_CNs(simulated_chromosomes):
    observed_CNs = {}
    for chrom_type in simulated_chromosomes:
        observed_CNs[chrom_type] = [len([x for x in simulated_chromosomes[chrom_type] if paternal == x["paternal"]]) for paternal in [True,False]]

    return(observed_CNs)

def count_CN_multiplicities(observed_CNs):
    multiplicities = {}
    
    for chrom_type in observed_CNs:
        for CN in observed_CNs[chrom_type]:
            if CN not in multiplicities: 
                multiplicities[CN] = 1

            else:
                multiplicities[CN] += 1

    return(multiplicities)

# for every copy number sum the precomputed values weighted against their multiplicity
# then adjust for CN’s not being able to go to zero
def CN_multiplicities_to_likelihoods(observed_CN_multiplicities):
    # these relative references will need to be modified before publishing 
    def CN_filename(CN):
        return precomputed_file_folder + "/collated_128_p128_v3_"+str(CN)+".npy"

    def ll(copy,multiplicity):
        return np.load(CN_filename(copy)) * multiplicity

    lls = None
    for copy in observed_CN_multiplicities:
        if lls is None:
            lls = ll(copy,observed_CN_multiplicities[copy])

        else:
            lls += ll(copy,observed_CN_multiplicities[copy])

    # account for the inability to lose all copies of a particular chromosome:
    likelihoods = np.exp(lls - np.log(1-np.exp(np.load(CN_filename(0)))) * observed_CN_multiplicities[0])

    # now we want to merge these computed likelihoods with the metadata columns:
    named_likelihoods = pkl.load(open(precomputed_file_folder+"/collated_128_p128_v3_list.pickle",'rb'))

    named_likelihoods.insert(3,"likelihood",likelihoods,True)
    named_likelihoods.columns = ["p_up","p_down","path","likelihood"]
    named_likelihoods.replace([np.inf,-np.inf], np.nan, inplace=True)
    # would like to investigate further as to what results in "np.inf"; there should be no likelihood of infinte value

    named_likelihoods.dropna(axis=0)
    named_likelihoods.sort_values(by=['likelihood'], inplace=True, ascending=False)

    total = np.nansum(named_likelihoods["likelihood"])
    print("total likelihood sum to normalise: "+str(total))
    named_likelihoods["likelihood"] /= total
    print("best likelihoods")
    print(named_likelihoods[:][0:300].to_string())

    return(named_likelihoods)


def likelihoods_to_marginal_likelihoods(likelihoods):
    marginal_likelihoods = likelihoods#.copy()
    marginal_likelihoods["mean_p_up"]   = marginal_likelihoods["p_up"]  *marginal_likelihoods["likelihood"]
    marginal_likelihoods["mean_p_down"] = marginal_likelihoods["p_down"]*marginal_likelihoods["likelihood"]
    marginal_likelihoods = marginal_likelihoods.groupby(["path"], as_index=False)[['likelihood','mean_p_up','mean_p_down']].sum()
    marginal_likelihoods.dropna(axis=0)
    marginal_likelihoods.sort_values(by=['likelihood'], inplace=True, ascending=False)
    total = sum(marginal_likelihoods["likelihood"]) 
    print("total marginal liklihoods sum to normalise: "+str(total))
    marginal_likelihoods["mean_p_up"] /= marginal_likelihoods["likelihood"] 
    marginal_likelihoods["mean_p_down"] /= marginal_likelihoods["likelihood"]
    print("best marginal likelihoods")
    print(marginal_likelihoods[:][0:20].to_string())
    print(sum(marginal_likelihoods["likelihood"]))

    return(marginal_likelihoods)


###### STEP 3; calculate the SNV multiplicities of each chromosome
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


##### STEP 4; now we have SNV counts, make all possible trees that could explain those SNV counts for the given epoch structure “(pre,mid, post)”
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
            new_CNs = (CN,tree[0]-CN)

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

    if len(tree) == 3:
        return((tree[0],complete_tree(tree[1]),complete_tree(tree[2])))

    elif len(tree) == 2:
        return((tree[0],complete_tree(tree[1]),complete_tree([tree[0]-tree[1][0]])))

    elif len(tree) == 1:
        if tree[0] == 0 or tree[0] ==1:
            return(tree)
        return((tree[0],complete_tree((tree[0]-int(tree[0]/2),)),complete_tree((int(tree[0]/2),))))
    else:
        assert(1==2) # throw a better error message than this
    return None


def complete_trees(trees):
    return([complete_tree(tree) for tree in trees])


# from the multiplicity counts of the chromosomes and the SNVs generate all the possible trees of every copy number:
def generate_trees(observed_CNs,SNV_CNs):
    SNV_CNs.sort(reverse = True)
    observed_CNs.sort(reverse = True)
    # print("SNV_CNs")
    # print(SNV_CNs)
    # print("observed_CNs")
    # print(observed_CNs)

    # initially we start with the following tree for each chromosome:
    trees = [(sum(observed_CNs),(max(observed_CNs),),(min(observed_CNs),))]

    for SNV_CN in SNV_CNs:
        # to save computational complexity we generate only trees from SNVs with copy numbers greater than 1
        if SNV_CN == 1:
            continue

        # insert the node at least once in every tree
        trees_with_new_node = insert_node(trees, SNV_CN)
 
        if trees_with_new_node == []:
            # then there isn’t anywhere left to insert the new node into these trees
            # this shouldn’t happen unless this SNV_CN is also in our observed_CNs
            assert(SNV_CN in observed_CNs)
            continue

        if SNV_CN in observed_CNs:
            # then it has already been inserted into the tree once in the first branch split
            # we still need to attempt inserting it again it might present in the alternative tree branch:
            trees = trees + trees_with_new_node
        else:
            # then we enforce that it is inserted at least once:
            trees = trees_with_new_node

        while(True):
            # continue to reinsert this node into more places in this tree until no new trees are found
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


##### STEP 5; now that all trees have been created, calculate all possible timings for each tree
##### 
##### 
##### 
##### 
##### 


def label_tree(tree, count, parents, label_to_copy):

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


def get_timings_per_tree(tree,epochs):
    # for a given tree assign labels to it and describe the parental relationships within the tree so that we can 
    # compute the timings of nodes in arrays

    labelled_tree, count, parents, label_to_copy = label_tree(copy.deepcopy(tree),0,{},{})
    # labelled_tree simply assigns a unique number to every node in the tree
    # count is the number of nodes in the tree
    # parents is a dictionary that describe the parental relationships between the unique numbers of the nodes in the tree
    # label_to_copy is a dictionary where the keys are the unique labels and the values are the copy number of the corresponding node

    count = count + 1

    # set up an empty array to begin with, each column will represent a particular node in the tree
    timings = np.array([None]*count)
    unique_tree_labels = range(count)

    for label in unique_tree_labels:
        if label == 0:
            # then it is the root node and there is no such thing as a timing for the first bifurcation
            # it is simply a natural split due to each person inheriting one of each CN
            timings = np.tile(timings, (1,1)) # this forces a potentially 1d array to be 2d
            timings[:,label] = 0 

        elif label_to_copy[str(label)] == 1:
            # then it is a leaf node and we can set it to be 
            timings = np.tile(timings, (1,1))
            timings[:,label] = epochs 
            # or label_to_copy[str(label)] == 0: 
            # I removed this second condition on forcing lost CNs to have thier BP prob calculated all the way to the end
            # it is a benefit to find out when they were most likely lost
            # it is also of potentially great computational cost so leave this comment here for further investigation in the future

        else:
            parent = int( parents[ str( label ) ] ) 

            # each row in timings is a particular solution to the problem of finding bifurcation times of each node in a given tree
            for row in range(len(timings)):
                parents_time = timings[row][parent]

                # not really too sure why parents_time could be none
                #assert(not (parents_time is None))
                if parents_time is None:
                    continue
                # don't know if this is the right thing, my quick look is that parents time will be none if the number of nestings are incompatible with the number of epochs.

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


def get_all_trees_and_timings(observed_SNV_multiplicities, observed_CNs):
    trees_and_timings = {}
    for chrom in observed_SNV_multiplicities:
        all_trees = generate_trees( 
            observed_CNs = observed_CNs[chrom],
            SNV_CNs = list(observed_SNV_multiplicities[chrom].keys())
                )

        epochs = pre*(pre>0) + mid*(mid>0) + post*(post>0) + (mid>=0) + (post>=0)

        trees_and_timings[chrom] = [get_timings_per_tree(x,epochs) for x in all_trees]
        trees_and_timings[chrom] = [x for x in trees_and_timings[chrom] if not None in x[3]]

        # this is potentially an error i haven't complettely investigated, 
        # it might be possible to have None in one row but not all rows of the array?

        if len(trees_and_timings[chrom]) == 0:
            print(trees_and_timings[chrom])
            #print("CAREFUL\n"*10)

        print(trees_and_timings[chrom])
    
    return(trees_and_timings)


##### STEP 6; now that all timing arrays for the nodes of each tree have been created, calculate the branch lengths
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
    CNs = [x for x in re.split("\(|\)|,|'", str(tree)) if x.isdigit()]
    unique_CNs = [int(x) for x in list(set(CNs))]
    unique_CNs = sorted(unique_CNs,reverse=True)
    unique_CNs = [str(x) for x in unique_CNs]

    for CN in unique_CNs:
        indices = find_indices(CNs,CN)
        new_stacked_branch_lengths = branch_lengths[:,indices].sum(axis=1)

        if CN == unique_CNs[0]:
            stacked_branch_lengths = new_stacked_branch_lengths

        else:
            stacked_branch_lengths = np.vstack((stacked_branch_lengths, new_stacked_branch_lengths))

    return((CNs, unique_CNs, branch_lengths, np.transpose(stacked_branch_lengths)))

# WHY DOES THIS NEED TO BE TRANSPOSED?

##### STEP 7; from the branch lengths calculate the BP likelihoods
##### 
##### 
##### 
##### 
##### 


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


def timing_struct_to_BP_likelihood_per_chrom(data, trees_and_timings, pre, mid, post):

    all_BP_likelihoods = []

    for these_tts in trees_and_timings:
        print(these_tts)
        if None in these_tts[3]:
            BP_likelihoods = -1

        else:
            CNs, unique_CNs, branch_lengths, stacked_branch_lengths = get_branch_lengths(these_tts)

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
            print(these_tts)
            print("copy numbers")
            print(CNs)
            print("branch lengths")
            print(branch_lengths)
            #>>> data = pickle.load(open("pre_mat129_u65_d10.precomputed.pickle",'rb'))
            # in the keys are all the possible paths...
            # each of these is a matrix that you can use to calculate the possible paths
        
            ends = these_tts[3]
            starts = ends - branch_lengths

            paths = np.zeros(ends.shape, dtype=object, order='C')
            likelihoods = np.zeros(ends.shape, dtype=float, order='C')

            for row in range(branch_lengths.shape[0]):
                for col in range(branch_lengths.shape[1]):
                    these_paths = path[starts[row][col]:ends[row][col]]
                    path_code = get_path_code(these_paths)

                    if CNs[col] == '1':
                        likelihood = data[path_code][1][1]

                    elif CNs[col] == '0':
                        likelihood = data[path_code][1][0]

                    else:
                        likelihood = data[path_code][1][2]

                    paths[row][col] = path_code
                    likelihoods[row][col] = likelihood

            print("starts")
            print(starts)
            print("ends")
            print(ends)
            print("paths")
            print(paths)
            print("likelihoods")
            print(likelihoods)

            ll = np.log(likelihoods)
            print("loglikelihoods")
            print(ll)
            BP_likelihoods = np.sum(ll[:,1:], axis=1)
            print("summed loglikelihoods")
            print(BP_likelihoods)

        all_BP_likelihoods += [BP_likelihoods]

    return(all_BP_likelihoods)


def get_BP_likelihoods(trees_and_timings,pre,mid,post,p_up,p_down):
    file = precomputed_file_folder + \
        "/precomputed/store_pre_pickle/pre_mat129_u"+str(int(p_up))+ \
        "_d"+str(int(p_down))+".precomputed.pickle"
    data = pkl.load(open(file,'rb'))
    BP_likelihoods = {}
    for chrom in trees_and_timings.keys():
        BP_likelihoods[chrom]  = timing_struct_to_BP_likelihood_per_chrom(
                data=data,
                trees_and_timings=trees_and_timings[chrom],
                pre=pre,
                mid=mid,
                post=post
                )
    return(BP_likelihoods)


##### STEP 8; from the branch lengths and the BP likelihoods calculate the join CN-SNV likelihoods
##### 
##### 
##### 
##### 
##### 


def get_poisson_loglikelihood(counts,stacked_branch_lengths,plambda,chrom,to_delete):
    A = np.log(stacked_branch_lengths.astype(float) * plambda * lengths[chrom]) * counts
    B = -np.tile( [scipy.special.gammaln(x+1) for x in counts], (stacked_branch_lengths.shape[0],1))
    C = -stacked_branch_lengths * plambda * lengths[chrom]

    summed = A + B + C
    summed = np.delete(summed,to_delete,1)
    total = np.sum(summed, axis=1)

    not_a_probability = any([x>0 for x in total])
    if not_a_probability:
        print(bits)
        print(A)
        print(B)
        print(C)
        print("\t"+str(summed))
        print("\t"+str(total))
        print(bits)
        print(plambda * lengths[chrom])
        assert(not not_a_probability)

    return(total)


def get_all_poisson_loglikelihoods_per_chr(timings,plambda,BP_likelihoods,observed_SNV_multiplicities,chrom): # these "timings" are on a per chromosome basis
    SNV_likelihoods = []
    for i in range(len(timings)):
        CNs, unique_CNs, branch_lengths, stacked_branch_lengths = get_branch_lengths(timings[i])

        counts = []
        for CN in unique_CNs:
            if int(CN) not in observed_SNV_multiplicities:
                counts += [0]
            else:
                counts += [observed_SNV_multiplicities[int(CN)]]


        if '0' not in unique_CNs:
            to_delete = [0]
        else:
            to_delete = [len(unique_CNs)-1]

        if unique_CNs[0] == '1' and unique_CNs[-1] == '0':
            to_delete = [0,len(unique_CNs)-1]

        #not_CN_0 = np.tile(np.array([int(x) > 0 for x in unique_CNs]),(len(BP_likelihoods[i]),1))
        #int_CNs = [int(x) for x in CNs]
        #if sum([max(int_CNs) == x for x in int_CNs]) == 1:
        #    not_root_CN = np.tile(np.array([int(x) != max(int_CNs) for x in unique_CNs]),(len(BP_likelihoods[i]),1))
        #    print(not_CN_0)
        #    not_CN_0 = not_root_CN * not_CN_0
        #    print(not_CN_0)

        #print(CNs)
        #print(unique_CNs)
        #print(branch_lengths)
        #print(stacked_branch_lengths)
        # put together BP and poisson likelihoods
        this_SNV_likelihood = get_poisson_loglikelihood(
                counts=counts, 
                stacked_branch_lengths=stacked_branch_lengths, 
                plambda=plambda,
                chrom=chrom,
                to_delete=to_delete
                ) 
        #print(chrom)
        #print(i)
        #print(this_SNV_likelihood)
        this_SNV_likelihood += BP_likelihoods[i]
        #print(this_SNV_likelihood)

        SNV_likelihoods += [this_SNV_likelihood]

    return(SNV_likelihoods)


def find_best_SNV_likelihood(plambda, timings, BP_likelihoods):
    SNV_likelihoods = {}
    best = {}

    chroms = timings.keys()

    for chrom in chroms:
        SNV_likelihoods[chrom] = get_all_poisson_loglikelihoods_per_chr(
                timings=timings[chrom], 
                plambda=plambda, 
                BP_likelihoods=BP_likelihoods[chrom],
                observed_SNV_multiplicities=observed_SNV_multiplicities[chrom],
                chrom=chrom # still need to pass in chrom for the lengths array
                )
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


def objective_function_SNV_loglik(plambda,timings,BP_likelihoods):
    total,best = find_best_SNV_likelihood(plambda,timings,BP_likelihoods)
    print(best)
    return(-total)


def path_code_to_pre_mid_post(path):
    bits = [int(x) for x in path.split("G")] + [-1,-1]
    pre, mid, post = bits[0:3]
    return((pre,mid,post))


##### STEP 9; write code to compare trees
##### 
##### 
##### 
##### 
##### 

# the following functions are mostly written with chatgpt3:

# a small function to see if two trees are identical:
# This function takes two trees as input, represented as lists, and returns a Boolean indicating whether they are topologically similar.
# The function first checks if the two trees have different lengths, in which case it returns False. 
# If the length is equal to 1, the function returns whether the single node value is equal. 
# If the length is greater than 1, the function recursively compares the two children trees.
def is_the_same_CN_tree(tree1,tree2):
    if len(tree1) != len(tree2):
        return False

    if len(tree1) == 1:
        return tree1[0] == tree2[0]

    return compare_trees(tree1[1], tree2[1]) and 
            compare_trees(tree1[2], tree2[2])


# a function that counts the number of matching nodes in the tree:
# The function now returns an integer instead of a Boolean, which represents the number of nodes in the two trees that are identical. 
# If the two trees have different lengths, the function returns 0. 
# If the length is equal to 1, the function returns 1 if the single node value is equal, or 0 otherwise. 
# If the length is greater than 1, the function recursively compares the two children trees and adds 1 if the current node value is equal.
def count_matching_CN_nodes(tree1,tree2):
    if len(tree1) != len(tree2):
        return 0
    if len(tree1) == 1:
        return 1 if tree1[0] == tree2[0] else 0
    return (compare_trees(tree1[1], tree2[1]) + compare_trees(tree1[2], tree2[2]) +
            (1 if tree1[0] == tree2[0] else 0))


# a function to count the number of nodes in just one tree:
# This function takes a single tree as input and returns the number of nodes in the tree. 
# If the length of the tree is equal to 1, the function returns 1, indicating that there's only one node in the tree. 
# If the length is greater than 1, the function recursively counts the number of nodes in the two children trees and adds 1 to represent the current node.
def count_nodes(tree):
    if len(tree) == 1:
        return 1

    return 1 + count_nodes(tree[1]) + count_nodes(tree[2])

# This function takes two trees as input, each represented as a dictionary, and returns a Boolean indicating whether they are topologically identical and have the same epoch_created value and copy_number value at each node. 
# The function first checks if the unique_identifier, epoch_created, and copy_number values are equal between the two trees. 
# If both trees have both child and complement keys missing, the function returns True. 
# If only one tree has a child or complement key missing, the function returns False. 
# If both trees have child or complement keys, the function recursively calls itself on the child or complement of the two trees. 
# If all checks return True, the function returns True, indicating that the two trees are topologically identical and have the same epoch_created value and copy_number value at each node.
def is_the_same_tree_by_epoch_and_time_created(tree1,tree2):
    #if tree1['unique_identifier'] != tree2['unique_identifier']:
    #    return False

    if tree1['epoch_created'] != tree2['epoch_created']:
        return False

    if tree1['copy_number'] != tree2['copy_number']:
        return False

    if tree1.get('child') is None and 
        tree2.get('child') is None and 
        tree1.get('complement') is None and 
        tree2.get('complement') is None:
        return True

    if (tree1.get('child') is None) != (tree2.get('child') is None):
        return False

    if (tree1.get('complement') is None) != (tree2.get('complement') is None):
        return False

    if tree1.get('child') is not None and not compare_trees(tree1['child'], tree2['child']):
        return False

    if tree1.get('complement') is not None and not compare_trees(tree1['complement'], tree2['complement']):
        return False

    return True


# Here's a modified version of the function that sums the absolute differences between the epoch_created values of each node for the nodes that have identical copy_number value:
# This function works similarly to the previous compare_trees function, but now adds the absolute difference between the epoch_created values of each node to sum if the copy_number values are equal. 
# The rest of the function remains the same as in the compare_trees function.
def compare_and_sum_trees(tree1, tree2):
    sum = 0

    #if tree1['unique_identifier'] != tree2['unique_identifier']:
    #    return 0

    if tree1['copy_number'] == tree2['copy_number']:
        sum += abs(tree1['epoch_created'] - tree2['epoch_created'])

    if tree1.get('child') is None and 
            tree2.get('child') is None and 
            tree1.get('complement') is None and 
            tree2.get('complement') is None:
        return sum

    if (tree1.get('child') is None) != (tree2.get('child') is None):
        return 0

    if (tree1.get('complement') is None) != (tree2.get('complement') is None):
        return 0

    if tree1.get('child') is not None:
        sum += compare_and_sum_trees(tree1['child'], tree2['child'])

    if tree1.get('complement') is not None:
        sum += compare_and_sum_trees(tree1['complement'], tree2['complement'])

    return sum


# write a function that takes a CN tree and a timings estimate and puts them together in the same structure as the truth data for comparison
# Here's a function that takes a tree data structure in the form [value, [child1, child2]] and converts it into a dictionary-like data structure in the form {'copy_number': value, 'child': child1, 'complement': child2}:
def convert_to_dict_tree(tree):
    copy_number, children = tree
    child1, child2 = children
    child1 = convert_to_dict_tree(child1) if isinstance(child1, list) else child1
    child2 = convert_to_dict_tree(child2) if isinstance(child2, list) else child2

    return {'copy_number': copy_number, 'child': child1, 'complement': child2}
# The function first unpacks the value and [child1, child2] components of the input tree and assigns them to the variables copy_number and children. It then further unpacks children into child1 and child2. If either child1 or child2 are lists, the function recursively calls convert_to_dict_tree


# now write a function with two arguments, the first is a dictionary-like tree data structure in the form {'copy_number': value, 'child': child1, 'complement': child2} called tree and the second is a list of numbers which is as long as the number of nodes in the tree. Add these numbers to the tree like data structure under the key "epoch_created" in a depth first way
# Here's a function that takes a dictionary-like tree data structure tree and a list of numbers and adds the numbers to the tree data structure under the key 'epoch_created' in a depth-first manner:
def add_epoch_created(tree, numbers):
    tree['epoch_created'] = numbers.pop(0)

    if 'child' in tree:
        add_epoch_created(tree['child'], numbers)

    if 'complement' in tree:
        add_epoch_created(tree['complement'], numbers)

    return tree

# The function first adds the first element of numbers to the tree dictionary under the key 'epoch_created' and removes it from numbers.
# If the 'child' key is present in the tree dictionary, the function recursively calls add_epoch_created on the value of the 'child' key with the updated numbers list. 
# If the 'complement' key is present, the function performs a similar operation. The final tree with the 'epoch_created' values is returned.

# now we write a function that can fully convert to the discovered tree and turn it into the form of the generated tree
def CN_tree_list_and_epoch_array_to_dictionary_tree(CN_tree,epoch_list):
    dict_tree = convert_to_dict_tree(CN_tree)
    dict_tree = add_epoch_created(dict_tree, epoch_list)

    return(dict_tree)



##### STEP 10; run the simulation 
#and try to find the parameters that created the simulation by optimising the likelihood of the simulated genome
##### 
##### 
##### 
##### 
##### 

print("START")
do_simulation = False #
do_simulation = True 
cache_results = False

pre = 2
mid = 2
post = -1
p_up=0.13
p_down=0.13
rate = 5

real_pre = pre
real_mid = mid
real_post = post
real_p_up = p_up
real_p_down = p_down
real_rate = rate


if do_simulation:
    simulated_chromosomes = simulate_single_with_poisson_timestamps_names(
            p_up=p_up, 
            p_down=p_down, 
            pre=pre, 
            mid=mid, 
            post=post, 
            rate=rate)

    truth_trees = create_truth_trees(simulated_chromosomes)
    for chrom_type in truth_trees:
        print(truth_trees[chrom_type])

    CN_trees = CN_trees_from_truth_trees(truth_trees)
    for chrom_type in CN_trees:
        print(CN_trees[chrom_type])

    # now we have to find a way to do the comparison of the trees
    # we want to see two things, 
    #   1) what percentage of the strucutre of the tree is correct, is it correct topologically?
    #   2) how correct are the timing estimates within the tree?
    #   3) what percentage of the relative timing estimates are correct?

    exit()

    print("observed chromosomal copynumbers")
    observed_CNs = count_CNs(simulated_chromosomes=simulated_chromosomes)
    print(observed_CNs)

    print("observed copynumber multiplicities")
    observed_CN_multiplicities = count_CN_multiplicities(observed_CNs=observed_CNs)
    print(observed_CN_multiplicities)

    print("likelihoods")
    likelihoods = CN_multiplicities_to_likelihoods(observed_CN_multiplicities=observed_CN_multiplicities)
    print(likelihoods)

    print("marginals")
    marginal_likelihoods = likelihoods_to_marginal_likelihoods(likelihoods=likelihoods)
    print(marginal_likelihoods)

    if cache_results:
        # need to save the important datastructures up to hear and then just work onwards from here to speed up development
        d = shelve.open('file.txt')           
        # in d is a dictionary type file that you can save variables:
        d['simulated_chromosomes'] = simulated_chromosomes
        d['observed_CNs'] = observed_CNs
        d['likelihoods'] = likelihoods
        d['marginal_likelihoods'] = marginal_likelihoods
        d.close()

if not do_simulation:
    # then load the most recently cached result:
    import shelve
    d = shelve.open('file.txt')
    simulated_chromosomes = d['simulated_chromosomes']
    observed_CNs = d['observed_CNs']
    likelihoods = d['likelihoods']
    marginal_likelihoods = d['marginal_likelihoods']
    d.close()


##### STEP 10; 
##### for each top result generate the potential trees and timing possibilities that explain that result. Use these trees and timings to come up with a better estimate of the likelihood
##### 
##### 
##### 
##### 



print("SNV multiplicities")
observed_SNV_multiplicities = count_SNV_multiplicities(simulated_chromosomes)
print(observed_SNV_multiplicities)

SEARCH_DEPTH = 1 
results = []
for res in range(SEARCH_DEPTH):
    path = marginal_likelihoods["path"].iloc[res]
    p_up = marginal_likelihoods['mean_p_up'].iloc[res].round()
    p_down = marginal_likelihoods['mean_p_down'].iloc[res].round()
    pre, mid, post = path_code_to_pre_mid_post(path)

    if res == SEARCH_DEPTH - 1:
        p_up = real_p_up*100
        p_down = real_p_down*100
        pre = real_pre
        mid = real_mid
        post = real_post

    print("calculate joint SNV-CN likelihood")
    print("pre: "+ str(pre))
    print("mid: "+str(mid))
    print("post: "+str(post))

    trees_and_timings = get_all_trees_and_timings(
            observed_SNV_multiplicities = observed_SNV_multiplicities,
            observed_CNs = observed_CNs
            )

    print("investigate the problem")
    for chrom in trees_and_timings:
        print(chrom)
        if len(trees_and_timings[chrom]) < 1:
            print(trees_and_timings[chrom])

    if res == SEARCH_DEPTH - 1:
        for chrom in trees_and_timings:
            assert(len(trees_and_timings[chrom]) >= 1)


    # need to take the rpinting of information out of the functions and into this main function only
    # if information needs to be printed then it needs to be able to be printed from this function here

    # from scipy.optimize import NonlinearConstraint, Bounds
    # output = scipy.optimize.minimize(fun=objective_function_SNV_loglik,x0=5,args=(timings))
    # output = scipy.optimize.differential_evolution(func=objective_function_SNV_loglik,bounds=Bounds([0.], [20.]),args=(timings))
    # need to reintroduce some sort of non array search for both the rate parameter AND p_up AND p_down

    # put branch lengths here, they are getting passed in to both comput ehte SNV likelihood and the BP likelihoods below
    # should not create the exact same data structure twice

    BP_likelihoods = get_BP_likelihoods(
            trees_and_timings=trees_and_timings,
            pre=pre,
            mid=mid,
            post=post,
            p_up=p_up,
            p_down=p_down
            )

    for chrom in BP_likelihoods:
        print("chrom: "+str(chrom))
        print("BP_L: "+str(BP_likelihoods[chrom]))

    # at some point evaluate the relative value of the likelihood contributed from the BP model to the likelihood contributed by the SNV model
    outputs = []
    for plambda in range(100):
        lam = plambda*10
        outputs += [(lam,objective_function_SNV_loglik(lam,trees_and_timings,BP_likelihoods))]

    #iterating though linearly doesn't seem to be too bad. It isn't fast but it isn't too bad, something to fix later
    # maybe we can add in the BP likelihoods now. 
    outputs.sort(key=lambda x: x[1])
    print("genome likelihoods vs lambda parameter: "+str(outputs))

    # now we need to modify it so that we are also printing out the best subtree for each tree
    result = find_best_SNV_likelihood(outputs[0][0],trees_and_timings,BP_likelihoods)
    print(result)
    results += [[result[0],pre,mid,post,p_up,p_down,result]]

for res in sorted(results):
    print(res)

        
