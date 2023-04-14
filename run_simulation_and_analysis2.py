# import the required libraries
import copy
import cProfile
import numpy as np
import pandas as pd
import pickle as pkl
import pstats
import re
import scipy
import signal
import shelve
import sys
from more_itertools import locate
from scipy.optimize import minimize_scalar

def pretty_print(x):
    max_len = 10000
    new_str = str(x)
    if len(new_str) > max_len:
        new_str = new_str[:max_len]
    print(new_str)
    return(None)


##### TIMEOUT LIMIT RUNTIME
# Define a function to handle timeouts
def timeout_handler(signum, frame):
    raise TimeoutError("The program took too long to run.")

# Set a timeout of 4 hours
timeout = 4*60*60

# Set the signal handler for SIGALRM (which is used for timeouts)
signal.signal(signal.SIGALRM, timeout_handler)
signal.alarm(timeout)

# check the health of this program, many potentially unused functions are being ddeveloped here
import ast
def find_unused_functions(source_code):
    # Parse the source code into an abstract syntax tree
    tree = ast.parse(source_code)

    # Create a set to store the names of all defined functions
    defined_functions = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            defined_functions.add(node.name)

    # Create a set to store the names of all called functions
    called_functions = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                called_functions.add(node.func.id)
            elif isinstance(node.func, ast.Attribute):
                called_functions.add(node.func.attr)

    # Return the set difference between defined and called functions
    return defined_functions - called_functions

with open("run_simulation_and_analysis2.py") as f:
    source_code = f.read()
unused_functions = find_unused_functions(source_code)
pretty_print("Unused functions:" +str(unused_functions))
# This function first parses the source code of the Python program into an abstract syntax tree (AST). 
# It then walks through the AST to identify all defined functions and all called functions. 
# Finally, it computes the set difference between the defined and called functions to find the functions that are defined but not used.



##### STEP 1; WRITE A FUNCTION THAT CAN SIMULATE A GENOME
#####
#####
#####
#####
#####

# precomputed files
precomputed_file_folder = "/vast/scratch/users/lmcintosh/GD2/GD/"

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
for chrom_type in lengths:
    lengths[chrom_type] *= 1/100


def count_paternity(chromosomes, paternal):
    """
    Count the number of chromosomes with a particular paternal type.
    """
    return len([chrom for chrom in chromosomes if chrom["paternal"] == paternal and not chrom["dead"]])

def check_all_chrs_are_unique(simulated_chromosomes):
    """
    Check if all chromosomes have unique identifiers.
    """
    ids = [chrom["unique_identifier"] for chrom_type in simulated_chromosomes.values() for chrom in chrom_type]
    return len(ids) == len(set(ids))

def check_expected_keys_in_simulated_chromosomes_present(simulated_chromosomes):
    """
    Check if all keys that are expected to be present in simulated chromosomes are actually present.
    """
    expected_keys = ["SNVs", "paternal", "epoch_created", "parent", "unique_identifier"]
    for chrom_type in simulated_chromosomes.values():
        for chrom in chrom_type:
            if not sorted(expected_keys) == sorted(chrom.keys()):
                return False
    return True

def get_ev_string(pre, mid, post):
    """
    Generate a list of characters representing the evolutionary path the cancer takes 
    """
    ev_str = []
    if pre > 0:
        ev_str += ["A"] * pre
    if mid > -1:
        ev_str += ["G"]
    if mid > 0:
        ev_str += ["A"] * mid
    if post > -1:
        ev_str += ["G"]
    if post > 0:
        ev_str += ["A"] * post
    return ev_str


#@profile
def simulate_single_with_poisson_timestamps_names(p_up,p_down,pre,mid,post,rate,agnostic=False):
    # this function simulates the evolution of a cancer genome with whole chromosome copy number changes and SNV's
    # pre, mid and post model the number of epochs / cell cycles that this genome undergoes between rounds of GD if they occur
    # pre is the number of epochs that occur before the first round of genome doubling, if there is a round of GD
    # mid is the number of epochs that occur after the first round of genome doubling and before the second round of genome doubling
    #   if there is no first round of genome doubling then mid is set to be -1 AND post is also -1
    # post is the number of epochs that occur after the second round of genome doubling
    #   if there is no second round of genome doubling then post is -1
    # rate is the base rate at which SNVs occur per epoch, 
    # this rate is scaled so that each chromosome on average will gain a number of SNVs proportional to its chromosomal length

    # SNV_count will count the number of SNVs and also provide a unique number to each new SNV created
    SNV_count = 0


    # set up the initial state of the genome
    simulated_chromosomes = {}
    for chrom_type in range(23):
        simulated_chromosomes[chrom_type] = [
                {
                    "unique_identifier" : chrom_type + x, 
                    # each chromosome gets a completely unique number to identify it in the simulation
                    "parent" : -1, 
                    # the value of parent is the the parent chromosomes unique identifier, 
                    # the root node of each chromosomal tree is specified to be -1
                    "epoch_created" : 0, 
                    # the epoch a chromosome is created, i
                    # both maternal and paternal chromosomes are created at the begining 
                    "paternal" : (x == 0), 
                    # True if paternal, False if maternal 
                    "SNVs":[], 
                    # a list of dicitonaries that describe the SNVs found on the chromosome, 
                    # a future extension might be to give these SNVs specific locaitons and simulate intra-chromosomal CN changes
                    "dead": False 
                    # if the chromosome is lost in the simulation we still need to record it to lineage tracing for the truth tree
                    } for x in (0,23)]

    # chrom_count will count the number of chromosomes created throughout the simulation and provide each with a unique identity
    chrom_count = 46

    # now simulate forward from the initial state

    # what is the total number of epochs?
    # turns out there is a complicated way to add these up and also an easy way:
    assert(pre*(pre>=0) + mid*(mid>=0) + post*(post>=0) + (mid>=0) + (post>=0) == pre+mid+post+2)

    # for each epoch in the total number of epochs:
    ev_sequence = get_ev_string(pre,mid,post)
    if len(ev_sequence) == 0:
        return(simulated_chromosomes)

    print("pre,mid,post:"+str((pre,mid,post)))
    print("ev_sequence:"+str(ev_sequence))

    for epoch,epoch_type in enumerate(ev_sequence):
        # SIMULATE SNVs
        # simulate and insert any new SNVs into the genome
        for chrom_type in simulated_chromosomes:
            for chrom in simulated_chromosomes[chrom_type]:
                if chrom["dead"]:
                    # then the chromosome has already been lost in the simulation => not appropriate to simulate SNVs on it
                    continue
                # generate a random number of SNVs proportional to the length of the genome:
                additional_SNV_count = np.random.poisson(rate * lengths[chrom_type],1)[0]

                # add these SNVs to the chromosome 
                for x in range(SNV_count + 1, SNV_count + additional_SNV_count + 1):
                    chrom["SNVs"] += [{"unique_identifier":str(x), "epoch_created":epoch}] #, "location":location}] 
                    # location can be added in later too 

                # update the SNV count
                SNV_count = SNV_count + additional_SNV_count 

        # if there is a genome doubling it has to be after post, so post cannot equal -1 if mid does
        # enforce this assertion as it may be easy to forget this later on:
        if mid == -1:
            assert(post == -1) 

        # SIMULATE GD
        # simulate the changes in copy number, but keep track of the SNVs
        # first see if this is a genome doubling round or an anueploidy round (they are mutually exclusive)
        if (mid != -1 and epoch == pre) or (post != -1 and epoch == pre+1+mid): 
            assert(epoch_type == "G")

            for chrom_type in simulated_chromosomes:
                new_chromosomes = []
                for chrom in simulated_chromosomes[chrom_type]:
                    if chrom["dead"]:
                        continue
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
            assert(epoch_type == "A")
            # this is a standard round of aneuploidy
            for chrom_type in simulated_chromosomes:
                if agnostic:
                    # generate simulated chromosomes without consistent probabilities

                    rnds = random.randint(0,3)
                    for rnd in range(rnds):
                        # randomly select a number of chromosomes to lose
                        chroms = simulated_chromosomes[chrom_type] 

                        for which in random.sample(range(len(chroms)), random.randint(0,len(chroms))):
                            chroms[which]["dead"] = True
            
                        # randomly select a number of chromosomes to gain after that 
                        # Add copies of a random number of items from the list back to itself
                        for which in random.choices(range(len(chrom)), k=random.randint(0,len(chroms))):
                            if not chroms[which]["dead"]:
                                new_chrom = copy.deepcopy(chroms[which])
                                new_chrom["epoch_created"] = epoch
                                new_chrom["parent"] = chroms[which]["unique_identifier"]
                                chrom_count += 1
                                new_chrom["unique_identifier"] = chrom_count
                                chroms.append(new_chrom)

                        # now do some GD:
                        # we only want 0 1 or 2 rounds of GD, so:
                        if rnd != rnds: 
                            new_chroms = []
                            for chrom in chroms: 
                                new_chrom = copy.deepcopy(chrom)
                                new_chrom["epoch_created"] = epoch
                                new_chrom["parent"] = chrom["unique_identifier"]
                                chrom_count += 1
                                new_chrom["unique_identifier"] = chrom_count
                                new_chroms.append(new_chrom)
                            chroms += new_chroms
                        
                        # update the list of chromosomes
                        simulated_chromosomes[chrom_type] = chroms


                else:
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
                            if chrom["dead"]:
                                new_chromosomes += [chrom]
                                continue

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
                                new_chromosome = copy.deepcopy(chrom)
                                new_chromosome["dead"] = True
                                new_chromosomes += [new_chromosome]

                        # ensure that there is always at least one copy of every chromosome
                        if len([x for x in new_chromosomes if not x["dead"]]) != 0:
                            break 
                   
                    # add those chromosomes into the genome
                    simulated_chromosomes[chrom_type] = new_chromosomes

    assert(pre+mid+post+2 == len(ev_sequence))
    # some quick final sanity checks:
    if post == 0 or (mid == 0 and post == -1):
        assert(ev_sequence[-1] == "G")
        # if a genome doubling round just occurred then the copy number of every chromosome will be even:
        for chrom_type in simulated_chromosomes: 
            if count_paternity(simulated_chromosomes[chrom_type],paternal=True) % 2 == 0:
                pretty_print(simulated_chromosomes[chrom_type])
                pretty_print(count_paternity(simulated_chromosomes[chrom_type],paternal=True))
                assert(count_paternity(simulated_chromosomes[chrom_type],paternal=True) % 2 == 0)
                
            if count_paternity(simulated_chromosomes[chrom_type],paternal=False) % 2 == 0:
                pretty_print(simulated_chromosomes[chrom_type])
                pretty_print(count_paternity(simulated_chromosomes[chrom_type],paternal=False))
                assert(count_paternity(simulated_chromosomes[chrom_type],paternal=False) % 2 == 0)

    for chrom_type in simulated_chromosomes:
        assert(len(simulated_chromosomes) != 0)

    # some basic checks about the simulated chromosomes:
    assert(check_all_chrs_are_unique(simulated_chromosomes))
    #assert(check_expected_keys_in_simulated_chromosomes_present(simulated_chromosomes))

    return(simulated_chromosomes)

##### STEP 1b; from the simulated genome create a tree
#####
#####
#####
#####
#####

# now make a structure to compare the truth tree to the found tree
#@profile
def insert_node_into_truth_tree(tree,node):
    #pretty_print("node")
    #pretty_print("\t"+str(node))
    #pretty_print("tree")
    #pretty_print("\t"+str(tree))
    assert(node["unique_identifier"] != tree["unique_identifier"])

    if node["parent"] == tree["unique_identifier"]:
        if tree["child"] == None:
            tree["complement"] = copy.deepcopy(tree)
            tree["complement"]["epoch_created"] = node["epoch_created"]

            # UNSURE, the following two lines i am unsure about
            if node["parent"] == -1: # the root node of every chromosomal tree has a unique identifier of 1 
                tree["complement"]["copy_number"] = 0
		# it seems that this copynumber get overwritten if another node is found to have a parent of -1 below

            # because tree is already a copy of itself if should already have child and complement set to None
            #assert(tree["complement"]["child"] == None)
            #assert(tree["complement"]["complement"] == None)

            tree["child"] = copy.deepcopy(node)
            tree["child"]["child"] = None
            tree["child"]["complement"] = None


        else:
            # then we need to insert into the complement if it is the correct node or further downt he tree
            if node["parent"] == -1: # then it is nested directly under the parent 
                assert(node["unique_identifier"] < 46)
                tree["complement"] = copy.deepcopy(node)
                tree["complement"]["child"] = None
                tree["complement"]["complement"] = None
                assert("copy_number" not in  tree["complement"])
            else:
                # this is what needs to be fixed i think?
                # if tree complement is
                if tree["complement"] is None:
                    tree["complement"] = copy.deepcopy(node)
                    tree["complement"]["child"] = None
                    tree["complement"]["complement"] = None
                else:
                    tree["complement"] = insert_node_into_truth_tree(tree["complement"],node)
    else:
        # insert it below child or complement (but not at that level)
        if tree["child"] != None:
            tree["child"] = insert_node_into_truth_tree(tree["child"],node)

        if tree["complement"] != None:
            tree["complement"] = insert_node_into_truth_tree(tree["complement"],node)

    return(tree)

def insert_node_into_truth_tree_v2(tree, node):
    assert(node["unique_identifier"] != tree["unique_identifier"])

    if node["parent"] == tree["unique_identifier"]:
        # then the node must be inserted here
        if tree["child"] is None:
            # either under child
            tree["complement"] = copy.deepcopy(tree)
            tree["complement"]["epoch_created"] = node["epoch_created"]

            tree["child"] = copy.deepcopy(node)
            tree["child"]["child"] = None
            tree["child"]["complement"] = None

        else:
            # or under the complement
            if node["parent"] == -1:
                # need to force it if maternal or paternal as the root node is ficticious
                assert(node["unique_identifier"] < 46)
                tree["complement"] = copy.deepcopy(node)
                tree["complement"]["child"] = None
                tree["complement"]["complement"] = None
                assert("copy_number" not in tree["complement"])

            else:
                # otherwise we insert it further down into the complement
                tree["complement"] = insert_node_into_truth_tree_v2(tree["complement"], node)

    else:
        # then we don't know where the node needs to be inserted and need to find that location:
        if tree["child"] is not None:
            tree["child"] = insert_node_into_truth_tree_v2(tree["child"], node)

        if tree["complement"] is not None:
            tree["complement"] = insert_node_into_truth_tree_v2(tree["complement"], node)

    return tree



# now that the truth tree is created for each chromosomes,
#   remove the SNV lists themselves, 
#   count the copynumber of each node and how many unique SNVs there are at that copy number.

# create a recursive function to insert the correct copy number at each node in the tree
#@profile
def add_copynumber_to_tree(tree):
    # handle the base case of losing one of the original paternal chromosomes
    if "copy_number" in tree:
        return(tree)

    if tree["child"] == None:
        assert(tree["complement"] == None)
        if tree["dead"]:
            tree["copy_number"] = 0
        else:
            tree["copy_number"] = 1

    else:
        tree["child"] = add_copynumber_to_tree(tree["child"])
        tree["complement"] = add_copynumber_to_tree(tree["complement"])
        tree["copy_number"] = tree["child"]["copy_number"] + tree["complement"]["copy_number"]

    return(tree)

# create a recursive function to find the correct number of SNVs at a particular copy at each node in the tree:
#@profile
def add_SNV_multiplicity_to_tree(tree):
    if tree["copy_number"] == 0:
        tree["SNV_multiplicity"] = None

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

#@profile
def remove_SNVs_from_tree(tree):
    tree.pop("SNVs")

    if tree["child"] != None:
        assert(tree["complement"] != None)
        tree["child"] = remove_SNVs_from_tree(tree["child"])
        tree["complement"] = remove_SNVs_from_tree(tree["complement"])

    return(tree)


def remove_dead_nodes(tree):
    if tree is None:
        return None

    if tree['child']:
        tree['child'] = remove_dead_nodes(tree['child'])
    if tree['complement']:
        tree['complement'] = remove_dead_nodes(tree['complement'])

    if tree['child'] is None and tree['complement'] is None and tree.get('dead', False) and tree['parent'] != -1:
        return None

    return tree

def are_all_descendants_zero(node):
    if node is None:
        return True

    if node["copy_number"] != 0:
        return False

    return are_all_descendants_zero(node["child"]) and are_all_descendants_zero(node["complement"])


def remove_redundant_parents(node):
    def are_all_descendants_zero(node):
        if node is None:
            return True

        if node["copy_number"] != 0:
            return False

        return are_all_descendants_zero(node["child"]) and are_all_descendants_zero(node["complement"])

    def check_and_remove_redundant_node(node, main_child, other_child):
        if main_child is not None and main_child["copy_number"] == node["copy_number"]:
            if not are_all_descendants_zero(other_child):
                raise ValueError("Invalid tree: child and complement copy_number do not add up to parent's copy_number")
            
            # Update epoch_created of the main_child to that of its parent
            main_child["epoch_created"] = node["epoch_created"]
            return main_child

        return node

    if node is None:
        return None

    # Recursively process child and complement nodes
    node["child"] = remove_redundant_parents(node["child"])
    node["complement"] = remove_redundant_parents(node["complement"])

    # Check for the conditions to remove the node, only if the node has a parent
    if node["parent"] is not None:
        node = check_and_remove_redundant_node(node, node["child"], node["complement"])
        node = check_and_remove_redundant_node(node, node["complement"], node["child"])

    return node


def with_epoch_killed(tree,max_epochs):
    """
    The function calculates the epoch_killed for each node in the given tree.

    :param tree: A dictionary representing a tree node with keys:
                 'child', 'complement', 'epoch_created', and 'epoch_killed'.
    :return: The modified tree with 'epoch_killed' values added.
    """

    if tree is None:
        return None

    child = tree.get('child')
    complement = tree.get('complement')

    if child is not None and complement is not None:
        child_epoch_created = child.get('epoch_created')
        complement_epoch_created = complement.get('epoch_created')

        if child_epoch_created != complement_epoch_created:
            raise ValueError("Epoch created values of child and complement do not match.")

        tree['child'] = with_epoch_killed(child,max_epochs)
        tree['complement'] = with_epoch_killed(complement,max_epochs)

        tree['epoch_killed'] = child['epoch_created']
    else:
        tree['epoch_killed'] = max_epochs 

    return tree


# The function takes in the tree as an argument and checks if the current node is dead or not. 
# If it is dead, the function moves the children of that node up to the parent of the dead node, by updating their parent field to be the same as the parent field of the dead node. 
# The function then recursively removes the dead node's children, returning either the child or complement of the dead node depending on which is not None. 
# If neither child is present, then the function returns None to remove the dead node from the tree. 
# If the current node is not dead, the function recursively removes the dead nodes from its children and returns the updated tree.
# The code checks if the parent of a dead node is equal to "-1" before removing it. If the parent is "-1", then the node is not removed and the function moves on to its children.

#@profile
def create_truth_trees(simulated_chromosomes,max_epochs):
    #there is a tree fro every chromosome
    trees = {}
    for chrom_type in simulated_chromosomes:
        # first sort the nodes by the order they need to be inserted in:
        sorted_list = sorted([(x["unique_identifier"],x) for x in simulated_chromosomes[chrom_type]])

        # create the root node of the tree for this chrom_type
        tree = {'unique_identifier':-1,
                'parent':None,
                'epoch_created':-1,
                'paternal':None,
                'child':None,
                'complement':None,
                'SNVs':[]
                } 

        # insert all nodes and add metadat to tree and simplify:
        for new_node in sorted_list:
            trees[chrom_type] = insert_node_into_truth_tree_v2(tree,new_node[1])  # UNSURE EDIT
        for i in range(5):
            pretty_print("##### chrom_type:"+str(chrom_type))
        pretty_print("with nodes inserted:")
        pretty_print(trees[chrom_type])
        trees[chrom_type] = add_copynumber_to_tree(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        pretty_print("with copynumber annotated:")
        pretty_print(trees[chrom_type])
        trees[chrom_type] = add_SNV_multiplicity_to_tree(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        pretty_print("with SNV multiplicity:")
        pretty_print(trees[chrom_type])
        trees[chrom_type] = remove_SNVs_from_tree(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        pretty_print("with SNVs removed from the tree:")
        pretty_print(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        trees[chrom_type] = remove_dead_nodes(trees[chrom_type])
        pretty_print("with dead nodes removed:")
        pretty_print(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        trees[chrom_type] = remove_redundant_parents(trees[chrom_type])
        pretty_print("with redundant nodes removed:")
        pretty_print(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        trees[chrom_type] = with_epoch_killed(trees[chrom_type],max_epochs)
        pretty_print("with epoch killed:")
        pretty_print(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        pretty_print(make_left_heavy(CN_tree_from_truth_tree(trees[chrom_type])))
        pretty_print("######\n"*5)

    return(trees)

#@profile
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
    if len(tree) == 1:
        return tree
    else:
        left_child = tree[1]
        right_child = tree[2] #if len(tree) == 3 else [0]

        if left_child[0] < right_child[0]:
            return [tree[0], make_left_heavy(right_child), make_left_heavy(left_child)]
        else:
            return [tree[0], make_left_heavy(left_child), make_left_heavy(right_child)]


#@profile
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
#@profile
def count_CNs(simulated_chromosomes):
    observed_CNs = {}
    for chrom_type in simulated_chromosomes:
        observed_CNs[chrom_type] = [len([x for x in simulated_chromosomes[chrom_type] if paternal == x["paternal"] and not x["dead"]]) for paternal in [True,False]]

    return(observed_CNs)

#@profile
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
#@profile
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
    if 0 not in observed_CN_multiplicities:
        likelihoods = np.exp(lls)
    else:
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
    #pretty_print("total likelihood sum to normalise: "+str(total))
    named_likelihoods["likelihood"] /= total
    #pretty_print("best likelihoods")
    #pretty_print(named_likelihoods[:][0:300].to_string())

    return(named_likelihoods)



#@profile
def likelihoods_to_marginal_likelihoods(likelihoods, top, default_paths):
    # Create a copy of the likelihoods DataFrame
    marginal_likelihoods = likelihoods.copy()

    # Calculate the mean probabilities
    marginal_likelihoods["mean_p_up"] = marginal_likelihoods["p_up"] * marginal_likelihoods["likelihood"]
    marginal_likelihoods["mean_p_down"] = marginal_likelihoods["p_down"] * marginal_likelihoods["likelihood"]

    # Group by path and sum likelihoods and mean probabilities
    marginal_likelihoods = marginal_likelihoods.groupby(["path"], as_index=False)[['likelihood','mean_p_up','mean_p_down']].sum()

    # Take the top n rows by likelihood value
    marginal_likelihoods = marginal_likelihoods.sort_values(by=['likelihood'], inplace=False, ascending=False)
    result = marginal_likelihoods.head(top)

    pretty_print("inside likelihoods to marginal likelihoods")
    pretty_print("marginal likelihoods")
    pretty_print(marginal_likelihoods)
    pretty_print("result")
    pretty_print(result)

    # Get the best row associated with each path in the default_paths list
    for path in default_paths:
        pretty_print("in marginals")
        pretty_print("path")
        pretty_print(path)
        best_row = marginal_likelihoods.loc[marginal_likelihoods["path"] == path]#.sort_values('likelihood', ascending=False).iloc[0]
        pretty_print("best_row")
        pretty_print(best_row)
        result = pd.concat([result,best_row])
        pretty_print("result")
        pretty_print(result)

    # Remove duplicate rows again, in case any of the default paths were already in the top n rows
    result = result.drop_duplicates(subset='path', keep='first')

    # Normalize the mean probabilities
    result["mean_p_up"] /= result["likelihood"]
    result["mean_p_down"] /= result["likelihood"]

    # Sort the result DataFrame by likelihood in descending order
    result = result.sort_values(by=['likelihood'], inplace=False, ascending=False)

    result["p_up"] = round(result["mean_p_up"])
    result["p_down"] = round(result["mean_p_down"])

    result = result[["p_up","p_down","path","likelihood"]]

    # Return the result DataFrame
    return result


#@profile
def likelihoods_to_best_likelihood_by_path(likelihoods, top, default_paths):
    # Create a copy of the likelihoods DataFrame
    df = likelihoods.copy()

    # Sort the DataFrame by likelihood in descending order
    df = df.sort_values('likelihood', ascending=False)

    # Remove duplicate rows
    df = df.drop_duplicates(subset='path', keep='first')

    # Take the top n rows by likelihood value
    result = df.head(top)
    for i in range(5):
        pretty_print('inside top likelihoods')
    # Get the best row associated with each path in the default_paths list
    for path in default_paths:
        pretty_print("in tops")
        pretty_print("path")
        pretty_print(path)
        pretty_print("result")
        pretty_print(result)
        pretty_print("best_row")
        path_rows = df.loc[df["path"] == path]
        best_row = path_rows.loc[path_rows['likelihood'].idxmax()].to_frame().T
        pretty_print("best_row")
        result = pd.concat([result,best_row])

    # Remove duplicate rows again, in case any of the default paths were already in the top n rows
    result = result.drop_duplicates(subset='path', keep='first')

    # Sort the result DataFrame by likelihood in descending order
    result = result.sort_values('likelihood', ascending=False)

    # Return the result DataFrame
    return result


###### STEP 3; calculate the SNV multiplicities of each chromosome
######
######
######
######

# first count the copy number of each SNV:
#@profile
def simulated_chromosomes_to_SNV_counts(simulated_chromosomes):
    SNV_copy_counter = {}
    for chrom_type in simulated_chromosomes:
        if chrom_type not in SNV_copy_counter:
            SNV_copy_counter[chrom_type] = {}

        for chrom in simulated_chromosomes[chrom_type]:
            if chrom["dead"]:
                continue
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

#@profile
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
#@profile
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
#@profile
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


#@profile
def complete_trees(trees):
    return([complete_tree(tree) for tree in trees])


# from the multiplicity counts of the chromosomes and the SNVs generate all the possible trees of every copy number:
#@profile
def generate_trees(observed_CNs,SNV_CNs):
    SNV_CNs.sort(reverse = True)
    observed_CNs.sort(reverse = True)
    pretty_print("SNV_CNs")
    pretty_print(SNV_CNs)
    pretty_print("observed_CNs")
    pretty_print(observed_CNs)

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
            pretty_print(SNV_CN)
            pretty_print(observed_CNs)
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


#@profile
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


#@profile
def get_timings_per_tree(tree,epochs):
    # for a given tree assign labels to it and describe the parental relationships within the tree so that we can 
    # compute the timings of nodes in arrays

    labelled_tree, count, parents, label_to_copy = label_tree(
            tree=copy.deepcopy(tree),
            count=0,
            parents={},
            label_to_copy={}
            )
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
            timings[:,label] =  0 # CHANGE BACK to -1 

###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this
###########WTF IS TIMINGS, need to understand this

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


#@profile
def get_all_trees_and_timings(observed_SNV_multiplicities, observed_CNs,pre,mid,post):
    trees_and_timings = {}
    for chrom in observed_SNV_multiplicities:
        pretty_print(chrom)

        all_trees = generate_trees( 
            observed_CNs = observed_CNs[chrom],
            SNV_CNs = list(observed_SNV_multiplicities[chrom].keys())
                )

        epochs = pre*(pre>0) + mid*(mid>0) + post*(post>0) + (mid>=0) + (post>=0)
        assert(epochs == pre+mid+post+2)

        trees_and_timings[chrom] = [get_timings_per_tree(x,epochs) for x in all_trees]
        trees_and_timings[chrom] = [x for x in trees_and_timings[chrom] if not None in x[3]]

        # this is potentially an error i haven't complettely investigated, 
        # it might be possible to have None in one row but not all rows of the array?

        if len(trees_and_timings[chrom]) == 0:
            pretty_print(trees_and_timings[chrom])
            #pretty_print("CAREFUL\n"*10)

        pretty_print(trees_and_timings[chrom])
    
    return(trees_and_timings)


##### STEP 6; now that all timing arrays for the nodes of each tree have been created, calculate the branch lengths
##### 
##### 
##### 
##### 
##### 

# this function will find the indices of the item_to_find in the list_to_check:
#@profile
def find_indices(list_to_check, item_to_find):
    indices = locate(list_to_check, lambda x: x == item_to_find)
    return list(indices)


# from the perspective of calculating the SNV likelihoods we care about how long it took for branches to bifurcate.
# we can obtain these branch_lengths from get_branch_lengths
#@profile
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
#@profile
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


#@profile
def timing_struct_to_BP_likelihood_per_chrom(data, trees_and_timings, pre, mid, post):

    all_BP_likelihoods = []

    for these_tts in trees_and_timings:
        pretty_print(these_tts)
        if None in these_tts[3]:
            BP_likelihoods = -1

        else:
            CNs, unique_CNs, branch_lengths, stacked_branch_lengths = get_branch_lengths(these_tts)

            path = []
            if pre > 0:
                path += ["A"]*pre
            if mid > -1:
                path += ["GD"]
            if mid > 0:
                path += ["A"]*mid
            if post > -1:
                path += ["GD"]
            if post > 0:
                path += ["A"]*post

            print("path:"+str(path))

            pretty_print("timings")
            pretty_print(these_tts)

            pretty_print("copy numbers")
            pretty_print(CNs)

            pretty_print("branch lengths")
            pretty_print(branch_lengths)

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
                    if path_code not in data:
                        print(path_code)
                        assert(path_code in data) # not sure why the path code is sometimes not in data...
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE
                    ######## CHECK BACK HERE

                    if CNs[col] == '1':
                        likelihood = data[path_code][1][1]

                    elif CNs[col] == '0':
                        likelihood = data[path_code][1][0]

                    else:
                        likelihood = data[path_code][1][2]

                    paths[row][col] = path_code
                    likelihoods[row][col] = likelihood

            pretty_print("starts")
            pretty_print(starts)
            pretty_print("ends")
            pretty_print(ends)
            pretty_print("paths")
            pretty_print(paths)
            pretty_print("likelihoods")
            pretty_print(likelihoods)

            ll = np.log(likelihoods)
            pretty_print("loglikelihoods")
            pretty_print(ll)
            BP_likelihoods = np.sum(ll[:,1:], axis=1)
            pretty_print("summed loglikelihoods")
            pretty_print(BP_likelihoods)

        all_BP_likelihoods += [BP_likelihoods]

    return(all_BP_likelihoods)


#@profile
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


#@profile
def get_poisson_loglikelihood(counts,stacked_branch_lengths,plambda,chrom,to_delete):
    A = np.log(stacked_branch_lengths.astype(float) * plambda * lengths[chrom]) * counts
    B = -np.tile( [scipy.special.gammaln(x+1) for x in counts], (stacked_branch_lengths.shape[0],1))
    C = -stacked_branch_lengths * plambda * lengths[chrom]

    summed = A + B + C
    summed = np.delete(summed,to_delete,1)
    total = np.sum(summed, axis=1)

    not_a_probability = any([x>0 for x in total])
    if not_a_probability:
        pretty_print(bits)
        pretty_print(A)
        pretty_print(B)
        pretty_print(C)
        pretty_print("\t"+str(summed))
        pretty_print("\t"+str(total))
        pretty_print(bits)
        pretty_print(plambda * lengths[chrom])
        assert(not not_a_probability)

    return(total)


#@profile
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
        #    pretty_print(not_CN_0)
        #    not_CN_0 = not_root_CN * not_CN_0
        #    pretty_print(not_CN_0)

        #pretty_print(CNs)
        #pretty_print(unique_CNs)
        #pretty_print(branch_lengths)
        #pretty_print(stacked_branch_lengths)
        # put together BP and poisson likelihoods
        this_SNV_likelihood = get_poisson_loglikelihood(
                counts=counts, 
                stacked_branch_lengths=stacked_branch_lengths, 
                plambda=plambda,
                chrom=chrom,
                to_delete=to_delete
                ) 
        #pretty_print(chrom)
        #pretty_print(i)
        #pretty_print(this_SNV_likelihood)
        this_SNV_likelihood += BP_likelihoods[i]
        #pretty_print(this_SNV_likelihood)

        SNV_likelihoods += [this_SNV_likelihood]

    return(SNV_likelihoods)


#@profile
def find_best_SNV_likelihood(plambda, timings, BP_likelihoods):
    print("plambda:"+str(plambda))
    print("timings:"+str(timings))
    print("BP_likelihoods:"+str(BP_likelihoods))
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

#@profile
def BP_and_SNV_loglik(plambda, p_up, p_down, trees_and_timings, pre, mid, post):
    
    BP_likelihoods = get_BP_likelihoods(
            trees_and_timings=trees_and_timings,
            pre=pre,
            mid=mid,
            post=post,
            p_up=p_up,
            p_down=p_down
            )

    total,best = find_best_SNV_likelihood(plambda,trees_and_timings,BP_likelihoods)

    return(-total)


#@profile
def find_BP_and_SNV_loglik(plambda_start, p_up_start, p_down_start, trees_and_timings, pre, mid, post, p_window, plambda_window):

    best_loglik = float("inf")
    best_p_up = 0
    best_p_down = 0
    best_plambda = 0

    #pretty_print(p_up_start)
    #pretty_print(p_down_start)
    #pretty_print(window)

    #for p_up in range(max(0,p_up_start - p_window), min(100,p_up_start + p_window + 1)):
    #    for p_down in range(max(0,p_down_start - p_window), min(100,p_down_start + p_window + 1)):
    # TEMP EDIT:
    for p_up in [p_up_start]:
        for p_down in [p_down_start]: 
            def optimize_func(plambda):
                return BP_and_SNV_loglik(plambda, p_up, p_down, trees_and_timings, pre, mid, post)

            #bounds = (0, None)
            bounds = (plambda_start*plambda_window, plambda_start/plambda_window)
            x0 = plambda_start


            res = minimize_scalar(optimize_func, bounds=bounds, options={'disp': True})
            if res.fun < best_loglik:
                best_loglik = res.fun
                best_p_up = p_up
                best_p_down = p_down
                best_plambda = res.x

                pretty_print(("Best Log-Likelihood:", best_loglik))
                pretty_print(("Best p_up:", best_p_up))
                pretty_print(("Best p_down:", best_p_down))
                pretty_print(("Best plambda:", best_plambda))

    # You can also use other optimization methods such as SLSQP, Nelder-Mead, Powell, COBYLA, TNC, BFGS, etc, as specified by the method argument.

    return(best_loglik, best_p_up, best_p_down, best_plambda, res)



#@profile
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
# now we have to find a way to do the comparison of the trees
#   1) what percentage of the strucutre of the tree is correct, is it correct topologically?
#   2) how correct are the timing estimates within the tree?
#   3) what percentage of the relative timing estimates are correct?


# the following functions in step 9 are simple typical tree comparison functions are were mostly written with the help of chatgpt3


# a small function to see if two trees are identical:
# This function takes two trees as input, represented as lists, and returns a Boolean indicating whether they are topologically similar.
# The function first checks if the two trees have different lengths, in which case it returns False. 
# If the length is equal to 1, the function returns whether the single node value is equal. 
# If the length is greater than 1, the function recursively compares the two children trees.
#@profile
def is_the_same_CN_tree(tree1,tree2):
    if len(tree1) != len(tree2):
        return False

    if len(tree1) == 1:
        return tree1[0] == tree2[0]

    return (is_the_same_CN_tree(tree1[1], tree2[1]) and 
            is_the_same_CN_tree(tree1[2], tree2[2]))


# a function to count the number of nodes in just one tree:
# This function takes a single tree as input and returns the number of nodes in the tree. 
# If the length of the tree is equal to 1, the function returns 1, indicating that there's only one node in the tree. 
# If the length is greater than 1, the function recursively counts the number of nodes in the two children trees and adds 1 to represent the current node.
#@profile
def count_nodes(tree):
    if "child" in tree and tree["child"] is not None:
        return count_nodes(tree["child"]) + count_nodes(tree["complement"]) + 1
    return 1 

# here's an example implementation of a function that takes a tree and swaps the 'child' and 'complement' nodes to ensure that 'child' always has a greater 'copy_number' than 'complement'. In case of ties, it also ensures that 'child' has a smaller 'epoch_created':
#@profile
def sort_tree_by_copy_number(tree):
    if tree is None:
        return None

    if 'child' in tree and tree['child'] is not None and \
        'complement' in tree and tree['complement'] is not None:
        if tree['complement'].get('copy_number') > tree['child'].get('copy_number'):
            tree['child'], tree['complement'] = tree['complement'], tree['child']
        elif tree['complement'].get('copy_number') == tree['child'].get('copy_number') and \
            tree['child'].get('epoch_index') > tree['complement'].get('epoch_index'):
            tree['child'], tree['complement'] = tree['complement'], tree['child']
        sort_tree_by_copy_number(tree['child'])
        sort_tree_by_copy_number(tree['complement'])
    elif 'child' in tree and tree['child'] is not None:
        sort_tree_by_copy_number(tree['child'])
    elif 'complement' in tree and tree['complement'] is not None:
        sort_tree_by_copy_number(tree['complement'])
    return tree

# This function uses recursion to traverse the tree and swaps the 'child' and 'complement' nodes whenever the 'copy_number' of 'complement' is greater than that of 'child', or when they have the same 'copy_number' and the 'epoch_created' of 'child' is greater than that of 'complement'. The function returns the modified tree.



# This function takes two trees as input, each represented as a dictionary, and returns a Boolean indicating whether they are topologically identical and have the same epoch_created value and copy_number value at each node. 
# The function first checks if the unique_identifier, epoch_created, and copy_number values are equal between the two trees. 
# If both trees have both child and complement keys missing, the function returns True. 
# If only one tree has a child or complement key missing, the function returns False. 
# If both trees have child or complement keys, the function recursively calls itself on the child or complement of the two trees. 
# If all checks return True, the function returns True, indicating that the two trees are topologically identical and have the same epoch_created value and copy_number value at each node.
#@profile
def is_the_same_dict_tree_by_epoch_and_time_created(tree1,tree2):
    if tree1['epoch_index'] != tree2['epoch_index']:
        return False

    if tree1['copy_number'] != tree2['copy_number']:
        return False

    if (tree1.get('child') is None and 
        tree2.get('child') is None and 
        tree1.get('complement') is None and 
        tree2.get('complement') is None):
        return True

    if (tree1.get('child') is None) != (tree2.get('child') is None):
        return False

    if (tree1.get('complement') is None) != (tree2.get('complement') is None):
        return False

    if (tree1.get('child') is not None and 
            not is_the_same_dict_tree_by_epoch_and_time_created(tree1['child'], tree2['child'])):
        return False

    if (tree1.get('complement') is not None and 
            not is_the_same_dict_tree_by_epoch_and_time_created(tree1['complement'], tree2['complement'])):
        return False

    return True


# Here's a modified version of the function that sums the absolute differences between the epoch_created values of each node for the nodes that have identical copy_number value:
# This function works similarly to the previous compare_trees function, but now adds the absolute difference between the epoch_created values of each node to sum if the copy_number values are equal. 
# The rest of the function remains the same as in the compare_trees function.
#@profile
def sum_tree_distance(tree1, tree2,diff_struct_is_inf=False):
    sum = 0

    if (tree1 is not None and tree2 is not None and
            tree1['copy_number'] == tree2['copy_number'] and 
            tree1["epoch_index"] is not None and 
            tree2["epoch_index"] is not None):
        sum += abs(tree1['epoch_index'] - tree2['epoch_index'])

    if tree1.get('child') is not None and tree2.get('child') is not None:
        sum += sum_tree_distance(tree1['child'], tree2['child'])
    elif diff_struct_is_inf:
        sum = float('inf')
    if tree1.get('complement') is not None and tree2.get('complement') is not None:
        sum += sum_tree_distance(tree1['complement'], tree2['complement'])
    elif diff_struct_is_inf:
        sum = float('inf')

    return sum

#@profile
def count_nodes_with_same_copy_number(tree1, tree2):
    count = 0
    if 'copy_number' in tree1 and 'copy_number' in tree2 and tree1['copy_number'] == tree2['copy_number']:
        count += 1
    for child_key in ['child', 'complement']:
        if child_key in tree1 and child_key in tree2 and tree1[child_key] is not None and tree2[child_key] is not None:
            count += count_nodes_with_same_copy_number(tree1[child_key], tree2[child_key])
    return count
#This function uses recursion to iterate through both trees and count the number of nodes that have the same 'copy_number' value. It first checks if the current nodes in both trees have a 'copy_number' key, and if they do and their values are equal, it increments the count by 1.
#The function then checks if both trees have child nodes with keys 'child' and 'complement', and if they do, it calls itself recursively for these child nodes and adds the returned count to the total count. 

#@profile
def count_nodes_with_same_properties(tree1, tree2):
    count = 0
    if 'copy_number' in tree1 and 'copy_number' in tree2 and tree1['copy_number'] == tree2['copy_number'] and \
    'epoch_index' in tree1 and 'epoch_index' in tree2 and tree1['epoch_index'] == tree2['epoch_index']:
        count += 1
    for child_key in ['child', 'complement']:
        if child_key in tree1 and child_key in tree2 and tree1[child_key] is not None and tree2[child_key] is not None:
            count += count_nodes_with_same_properties(tree1[child_key], tree2[child_key])
    return count
# This function is very similar to the previous one, but now it also checks if both 'copy_number' and 'epoch_created' fields are the same for the nodes in both trees. The function uses the same recursion strategy to traverse both trees.

# write a function that takes a CN tree and a timings estimate and puts them together in the same structure as the truth data for comparison
# Here's a function that takes a tree data structure in the form [value, [child1, child2]] and converts it into a dictionary-like data structure in the form {'copy_number': value, 'child': child1, 'complement': child2}:
#@profile
def convert_to_dict_tree(tree):
    if len(tree) == 3:
        copy_number, child1, child2 = tree
        child1 = convert_to_dict_tree(child1) if isinstance(child1, tuple) else child1
        child2 = convert_to_dict_tree(child2) if isinstance(child2, tuple) else child2
        return {'copy_number': copy_number, 'child': child1, 'complement': child2}
    elif len(tree) == 1:
        copy_number = tree[0]
        return {'copy_number': copy_number}

    else:
        sys.exit()
# The function first unpacks the value and [child1, child2] components of the input tree and assigns them to the variables copy_number and children. It then further unpacks children into child1 and child2. If either child1 or child2 are lists, the function recursively calls convert_to_dict_tree

# this function converts a tree in dictionary format to a list representation so its easy to visually inspect
#@profile
def convert_dict_tree_to_list(tree,total_epochs=None,is_truth=False):
    if tree is None:
        return None

    copy_number = tree.get('copy_number')
    epoch_index = tree.get('epoch_index')

    if is_truth:
        if 'child' in tree and tree["child"] is not None:
            epoch_index = tree["child"]["epoch_index"]
        else:
            epoch_index = total_epochs 

    child_tree = convert_dict_tree_to_list(tree.get('child'),total_epochs,is_truth)
    complement_tree = convert_dict_tree_to_list(tree.get('complement'),total_epochs,is_truth)
    if child_tree is None:
        return [(copy_number, epoch_index)]
    return [(copy_number, epoch_index), child_tree, complement_tree]

# This function uses recursion to traverse the dictionary tree and convert each node into a tuple containing the 'copy_number' and 'epoch_created' fields, and two lists representing the 'child' and 'complement' subtrees, respectively. If a node does not have a 'child' or 'complement' subtree, the corresponding list will be None.



# now write a function with two arguments, the first is a dictionary-like tree data structure in the form {'copy_number': value, 'child': child1, 'complement': child2} called tree and the second is a list of numbers which is as long as the number of nodes in the tree. Add these numbers to the tree like data structure under the key "epoch_created" in a depth first way
# Here's a function that takes a dictionary-like tree data structure tree and a list of numbers and adds the numbers to the tree data structure under the key 'epoch_created' in a depth-first manner:
#@profile
def add_epoch_bifurcated(dict_tree, epoch_list):
    #pretty_print(dict_tree)
    #pretty_print(epoch_list)
    dict_tree['epoch_created'] = epoch_list.pop(0)

    if 'child' in dict_tree:
        add_epoch_bifurcated(dict_tree['child'], epoch_list)

    if 'complement' in dict_tree:
        add_epoch_bifurcated(dict_tree['complement'], epoch_list)

    return dict_tree

# The function first adds the first element of numbers to the tree dictionary under the key 'epoch_created' and removes it from numbers.
# If the 'child' key is present in the tree dictionary, the function recursively calls add_epoch_created on the value of the 'child' key with the updated numbers list. 
# If the 'complement' key is present, the function performs a similar operation. The final tree with the 'epoch_created' values is returned.


def update_to_epoch_created(tree, parent_epoch=None):
    if 'child' in tree and tree['child'] is not None:
        update_to_epoch_created(tree['child'], tree['epoch_created'])
    
    if 'complement' in tree and tree['complement'] is not None:
        update_to_epoch_created(tree['complement'], tree['epoch_created'])

    if parent_epoch is not None:
        tree['epoch_created'] = parent_epoch

    return(tree)


# now we write a function that can fully convert to the discovered tree and turn it into the form of the generated tree
#@profile
def CN_tree_list_and_epoch_array_to_dictionary_tree(CN_tree,epoch_list):
    dict_tree = convert_to_dict_tree(CN_tree)
    #pretty_print(dict_tree)
    epoch_list = list(epoch_list)
    dict_tree = add_epoch_bifurcated(dict_tree, epoch_list)
    dict_tree = update_to_epoch_created(dict_tree)
    return dict_tree

def filter_tree(tree, keys_to_keep):
    filtered_tree = {}
    for key in keys_to_keep:
        if key in tree:
            filtered_tree[key] = tree[key]

    if 'child' in tree and tree['child'] is not None:
        filtered_tree['child'] = filter_tree(tree['child'], keys_to_keep)

    if 'complement' in tree and tree['complement'] is not None:
        filtered_tree['complement'] = filter_tree(tree['complement'], keys_to_keep)

    return filtered_tree

##### STEP 10; run the simulation 
#and try to find the parameters that created the simulation by optimising the likelihood of the simulated genome
##### 
##### 
##### 
##### 
##### 

pretty_print("START")
test_case = sys.argv[1]
#test_case = 10000 #sys.argv[1]
do_simulation = False 
do_simulation = True 

cache_results = False
cache_results = True

do_search = False
do_search = True

import random
import math

def random_decimal(min_value, max_value, decimal_places):
    random_number = random.uniform(min_value, max_value)
    return round(random_number, decimal_places)

def random_integer_log_scale(min_value, max_value):
    log_min = math.log(min_value)
    log_max = math.log(max_value)
    random_log = random.uniform(log_min, log_max)
    return int(round(math.exp(random_log)))

def biased_sample(p,min_value,max_value):
    other_values_count = max_value - min_value
    
    if random.random() < p:
        return min_value
    else:
        return random.randint(min_value + 1, max_value)


# the parameters that govern the search depth:
top = 2# top describes how many of the top solutions to go through
p_window = 0
plambda_window = 0.1

pretty_print("SEARCH PARAMETERS ARE: ")
pretty_print("BP search depth: "+ str(top))
pretty_print("additive window width to search around top estimates of enuploidy probabilities: " + str(p_window))
pretty_print("multiplicative window width to search around top estimates of the poisson parameter: " + str(p_window))


max_default_path_length = 0 #2
default_paths = [str(x) for x in range(max_default_path_length)]
default_paths += [str(x) + "G" + str(y) 
        for x in range(max_default_path_length) 
        for y in range(max_default_path_length) 
        if x+y <= max_default_path_length]

pretty_print("default paths that will always get searched through")
pretty_print(default_paths)
pretty_print("default path lengths")
pretty_print(len(default_paths))


if do_simulation:
    pretty_print("Doing simulation")
    # the parameters that generated the simulation:
    max_epochs = 4
    pre = random.randint(0,max_epochs) 
    mid = biased_sample(0.5,-1,max_epochs) #-1 
    post = biased_sample(0.8,-1,max_epochs)

    if mid == -1 and post >=0:
        temp = mid
        mid = -1
        post = temp


    total_epochs = pre+mid+post+(pre>=0)+(mid>=0)
    p_up = random_decimal(0.10,0.30,2) #0.13
    p_down = random_decimal(p_up*0.5,min(0.3,p_up*2),2) #0.13
    rate = 10 #random_integer_log_scale(100,100000) #100
    # a low rate helps for debugging purposes...

    pretty_print("SIMULATION PARAMETERS ARE: ")
    pretty_print("pre: "+ str(pre))
    pretty_print("mid: "+ str(mid))
    pretty_print("post: " + str(post))
    pretty_print("p_up: " + str(p_up))
    pretty_print("p_down: " + str(p_down))
    pretty_print("rate: " + str(rate))
    pretty_print("")


    # save the true parameters because there are name collisions:
    real_pre = pre
    real_mid = mid
    real_post = post
    real_p_up = p_up
    real_p_down = p_down
    real_rate = rate


    simulated_chromosomes = simulate_single_with_poisson_timestamps_names(
            p_up=p_up, 
            p_down=p_down, 
            pre=pre, 
            mid=mid, 
            post=post, 
            rate=rate)

    pretty_print("Simulated genome was:")
    for chrom in simulated_chromosomes:
        pretty_print("chrom: " + str(chrom))
        pretty_print(simulated_chromosomes[chrom])

    pretty_print("the truth trees are:")
    truth_trees =  create_truth_trees(simulated_chromosomes=simulated_chromosomes,max_epochs=pre+mid+post+2)
    for chrom_type in truth_trees:
        pretty_print(truth_trees[chrom_type])

    pretty_print("the CN simplified trees are:")
    CN_trees = CN_trees_from_truth_trees(truth_trees)
    for chrom_type in CN_trees:
        pretty_print(CN_trees[chrom_type])


    pretty_print("observed chromosomal copynumbers")
    observed_CNs = count_CNs(simulated_chromosomes=simulated_chromosomes)
    pretty_print(observed_CNs)

    pretty_print("observed copynumber multiplicities")
    observed_CN_multiplicities = count_CN_multiplicities(observed_CNs=observed_CNs)
    pretty_print(observed_CN_multiplicities)

    pretty_print("likelihoods")
    likelihoods = CN_multiplicities_to_likelihoods(observed_CN_multiplicities=observed_CN_multiplicities)
    pretty_print(likelihoods)

    pretty_print("marginals")
    marginal_likelihoods = likelihoods_to_marginal_likelihoods(
            likelihoods = likelihoods,
            top = top,
            default_paths = default_paths)

    pretty_print(marginal_likelihoods)

    pretty_print("top likelihoods")
    top_likelihoods = likelihoods_to_best_likelihood_by_path(
            likelihoods = likelihoods,
            top = top,
            default_paths = default_paths)

    pretty_print(top_likelihoods)


    marginal_likelihoods = marginal_likelihoods.reindex(columns=top_likelihoods.columns)
    searchable_likelihoods = pd.concat([top_likelihoods,marginal_likelihoods], ignore_index=True)
    searchable_likelihoods = searchable_likelihoods.drop_duplicates(subset=['path','p_up','p_down'], keep='first')
    searchable_likelihoods = searchable_likelihoods.sort_values(by='path')
    #searchable_likelihoods = searchable_likelihoods[searchable_likelihoods['likelihood'] > 1e-9]

    pretty_print("The best likelihoods are:")
    pretty_print(searchable_likelihoods)
    pretty_print("paths to search: "+str(len(searchable_likelihoods.index)))


    if cache_results:
        # need to save the important datastructures up to here and then just work onwards from here to speed up development
        d = shelve.open('file_'+str(test_case)+'.txt')           
        # in d is a dictionary type file that you can save variables:
        d['simulated_chromosomes'] = simulated_chromosomes
        d['truth_trees'] = truth_trees
        d['CN_trees'] = CN_trees
        d['observed_CNs'] = observed_CNs
        d['observed_CN_multiplicities'] = observed_CN_multiplicities
        d['likelihoods'] = likelihoods
        d['marginal_likelihoods'] = marginal_likelihoods
        d['top_likelihoods'] = top_likelihoods
        d['searchable_likelihoods'] = searchable_likelihoods
        d['pre'] = pre
        d['mid'] = mid
        d['post'] = post
        d['p_up'] = p_up
        d['p_down'] = p_down
        d['rate'] = rate
        d.close()

if not do_simulation:
    # then load the most recently cached result:
    #d = shelve.open('file.txt')
    d = shelve.open('file_'+str(test_case)+'.txt')           
    simulated_chromosomes = d['simulated_chromosomes']
    truth_trees = d['truth_trees']
    CN_trees = d['CN_trees']
    observed_CNs = d['observed_CNs']
    observed_CN_multiplicities = d['observed_CN_multiplicities']
    likelihoods = d['likelihoods']
    marginal_likelihoods = d['marginal_likelihoods']
    top_likelihoods = d['top_likelihoods']
    searchable_likelihoods = d['searchable_likelihoods']
    pre = d['pre']
    real_pre = pre
    mid = d['mid']
    real_mid = mid
    post = d['post']
    real_post = post
    p_up = d['p_up']
    real_p_up = p_up
    p_down = d['p_down']
    real_p_down = p_down
    rate = d['rate']
    real_rate = rate


    pretty_print("simulated_chromsomes")
    pretty_print(simulated_chromosomes)

    pretty_print("truth_trees")
    pretty_print(truth_trees)

    pretty_print("CN_trees")
    pretty_print(CN_trees)

    pretty_print("observed_CNs")
    pretty_print(observed_CNs)

    pretty_print("observed_CN_multiplicities")
    pretty_print(observed_CN_multiplicities)
    
    pretty_print("likelihoods")
    pretty_print(likelihoods)
    
    pretty_print("marginal_likelihoods")
    pretty_print(marginal_likelihoods)

    pretty_print("top_likelihoods")
    pretty_print(top_likelihoods)
    
    pretty_print("searchable_likelihoods")
    pretty_print(searchable_likelihoods)
    d.close()


##### STEP 10; 
##### for each top result generate the potential trees and timing possibilities that explain that result. Use these trees and timings to come up with a better estimate of the likelihood
##### 
##### 
##### 
##### 

#@profile
def sum_SNV_counts(observed_SNV_multiplicities):
    d = observed_SNV_multiplicities
    total = 0
    for key1 in d:
        for key2 in d[key1]:
            total += d[key1][key2]
    return total


#@profile
def sum_chrom_multiplicities(observed_CN_multiplicities):
    return sum(observed_CN_multiplicities.values())



pretty_print("SNV multiplicities")
observed_SNV_multiplicities = count_SNV_multiplicities(simulated_chromosomes)
pretty_print(observed_SNV_multiplicities)

if do_search:
    pretty_print("searchable likelihoods")
    pretty_print(searchable_likelihoods)

    SEARCH_DEPTH = len(searchable_likelihoods)
    #SEARCH_DEPTH = 0
    results = []
    all_trees_and_timings = []
    for res in range(SEARCH_DEPTH+1):
        print("############################\n"*10)
        print(res)

        if res == SEARCH_DEPTH:
            p_up = int(real_p_up*100)
            p_down = int(real_p_down*100)
            pre = real_pre
            mid = real_mid
            post = real_post
            path = get_ev_string(pre,mid,post) 

        else:
            path = searchable_likelihoods["path"].iloc[res]
            p_up = int(searchable_likelihoods['p_up'].iloc[res])
            p_down = int(searchable_likelihoods['p_down'].iloc[res])
            pre, mid, post = path_code_to_pre_mid_post(path)

        pretty_print("path: "+str(path))
        pretty_print("p_up: "+str(p_up))
        pretty_print("p_down: "+str(p_down))
        pretty_print("pre: "+str(pre))
        pretty_print("mid: "+str(mid))
        pretty_print("post: "+str(post))

        trees_and_timings = get_all_trees_and_timings(
            observed_SNV_multiplicities = observed_SNV_multiplicities,
            observed_CNs = observed_CNs,
            pre=pre,
            mid=mid,
            post=post
            )

        pretty_print("investigate the problem")
        print("print the trees and timings info")
        for chrom in trees_and_timings:
            pretty_print(chrom)
            pretty_print(trees_and_timings[chrom])

        #for chrom in trees_and_timings:
        #    assert(len(trees_and_timings[chrom]) >= 1)
        # surprisingly this may not be true?
        for chrom in trees_and_timings:
            if len(trees_and_timings[chrom]) == 0:
                continue



        pretty_print("SELECT the best estimates")
        #total_SNV_time = 0 # rtotal time refers to the number of epochs under which SNVs evolve.
        #if pre > 0:
        #    total_SNV_time += pre
        #if mid > 0:
        #    total_SNV_time += mid
        #if post > 0:
        #    total_SNV_time += post
        total_SNV_time = pre+mid+post+2+1 # SNVs occur every round not just in anueploidy rounds so the above commented code is wrong
        # BUT it also occurs in the first round when nothing happens, so we also add a plus one

        total_SNVs = sum_SNV_counts(observed_SNV_multiplicities)
        pretty_print("total SNVs")
        pretty_print(total_SNVs)
        pretty_print("observed_SNV_multiplicities")
        pretty_print(observed_SNV_multiplicities)

        total_chromosomes = sum_chrom_multiplicities(observed_CN_multiplicities)
        plambda_start = float(total_SNVs) / float(total_SNV_time) / float(total_chromosomes) * 23
        pretty_print("plambda_start: " + str(plambda_start))  
        # this is a very rough estimate of expected value but it should align with the optimised for value. 
        # and this should be checked at some point and this line removed once checked

        # these need to be changed to the values learnt in the array.
        p_up_start = p_up
        p_down_start = p_down

        BP_SNV_output = find_BP_and_SNV_loglik(
                plambda_start = plambda_start, 
                p_up_start = p_up_start, 
                p_down_start = p_down_start, 
                trees_and_timings = trees_and_timings, 
                pre = pre, 
                mid = mid, 
                post = post, 
                p_window =p_window,
                plambda_window = plambda_window
                )

        best_loglik, best_p_up, best_p_down, best_plambda, result = BP_SNV_output
        print("BP_SNV_output")
        print(BP_SNV_output)
        
        results += [[best_loglik, pre, mid, post, best_p_up, best_p_down, best_plambda, result]]
        all_trees_and_timings += [trees_and_timings] 
        #d = shelve.open('file2.txt')
        d = shelve.open('file2_'+str(test_case)+'.txt')           
        d['results'] = results
        d['trees_and_timings'] = all_trees_and_timings
        d.close()
else:
    # load the results array
    #d = shelve.open('file2.txt')
    d = shelve.open('file2_'+str(test_case)+'.txt')           
    results = d['results']
    all_trees_and_timings = d['trees_and_timings']
    d.close()
    # at some point evaluate the relative value of the likelihood contributed from the BP model to the likelihood contributed by the SNV model

print(results)
print(all_trees_and_timings)


def order_tree_keys_alphabetically(tree):
    ordered_tree = {}
    for key in sorted(tree.keys()):
        if key in ['child', 'complement']:
            ordered_tree[key] = order_tree_keys_alphabetically(tree[key]) if tree[key] else None
        else:
            ordered_tree[key] = tree[key]
    return ordered_tree

def create_epoch_index(tree, key_from, key_to):
    """
    The function replaces the key called 'key_from' in the dictionary 'tree' with a new key called 'key_to',
    keeping the original value.

    :param tree: A dictionary representing a tree node.
    :param key_from: The original key to be replaced.
    :param key_to: The new key that will replace the original key.
    :return: The modified tree with the key replaced.
    """

    if key_from in tree:
        tree[key_to] = tree[key_from]
        del tree[key_from]

    if 'child' in tree and tree['child'] is not None:
        tree['child'] = create_epoch_index(tree['child'], key_from, key_to)

    if 'complement' in tree and tree['complement'] is not None:
        tree['complement'] = create_epoch_index(tree['complement'], key_from, key_to)

    return tree


all_results = {}
for result_index,res in enumerate(sorted(results)):
    pretty_print(res)
    val, pre_est, mid_est, post_est, p_up_est, p_down_est, plambda_est, result = res
    print(res)
    trees_and_timings = all_trees_and_timings[result_index]
    print(trees_and_timings)
    BP_likelihoods = get_BP_likelihoods(
            trees_and_timings = trees_and_timings,
            pre = pre_est,
            mid = mid_est,
            post = post_est,
            p_up = p_up_est,
            p_down = p_down_est
            )

    total, best = find_best_SNV_likelihood(plambda_est,trees_and_timings,BP_likelihoods)

    simulated_trees = create_truth_trees(simulated_chromosomes=simulated_chromosomes,max_epochs=pre_est+mid_est+post_est+2)
    simplified_simulated_trees = {}
    estimated_trees = {}
    for chrom in best:
        pretty_print("printing the best likelihoods found in the search")
        for i in range(10):
            pretty_print("chrom: "+str(chrom))
        max_lik, tree_index, row_index = best[chrom]

        pretty_print("trees and timings")
        pretty_print(trees_and_timings[chrom])
        pretty_print("best one")
        pretty_print(tree_index)
        pretty_print(trees_and_timings[chrom][tree_index])
        CN_tree, labelled_tree, count, timings, parents = trees_and_timings[chrom][tree_index]

        pretty_print("CN tree")
        pretty_print(CN_tree)

        pretty_print("labelled tree")
        pretty_print(labelled_tree)

        pretty_print("count")
        pretty_print(count)

        pretty_print("timings")
        pretty_print(timings)

        pretty_print("row_index")
        pretty_print(row_index)

        pretty_print("epoch_list")
        epoch_list = timings[row_index]
        pretty_print(epoch_list)

        pretty_print("parents")
        pretty_print(parents)

        estimated_tree = CN_tree_list_and_epoch_array_to_dictionary_tree(CN_tree,epoch_list)
        estimated_trees[chrom] = estimated_tree

        # epoch_created in each tree in 'estimated_trees' is actually the bifurcation times to create the next node. 
        # to be able to accurately compare the timing structre of estimated and true cancer trees we need to change one of these
        # epoch_created in the truth trees is actually the real time that the chromosome was created. 
        # In the estimated trees, the time the parent bifurcates is the time that the tree was actually created. 

        pretty_print("estimated tree")
        estimated_trees[chrom] = order_tree_keys_alphabetically(estimated_trees[chrom])
        pretty_print(estimated_trees[chrom])

        pretty_print("simulated tree")
        simplified_simulated_trees[chrom] = order_tree_keys_alphabetically(filter_tree(simulated_trees[chrom],keys_to_keep = ["copy_number","epoch_killed"]))
        pretty_print(simplified_simulated_trees[chrom])
        
    # now compare the results to the truth!
    # CN tree is a list that is estimated when optimising the lieklihood. 
    # Need a way to pass that back
    # epoch list is the row in the epoch array that was optimised for
    # also need a way to pass that back to here.
    # also need a way to annotate the estimated tree with estimated number of SNVs under each branch of the tree.
    # being able to estimate the individual timing of each SNV will be very valuable

    simulated_trees = simplified_simulated_trees #truth_trees

    total_nodes = 0
    num_chrom_with_correct_CN = 0
    num_chrom_with_correct_CN_and_epoch_created = 0
    average_distance_from_truth_of_epoch_created = 0
    for chrom in estimated_trees:
        simulated_trees[chrom] = create_epoch_index(tree=simulated_trees[chrom],key_from="epoch_killed",key_to="epoch_index")
        estimated_trees[chrom] = create_epoch_index(tree=estimated_trees[chrom],key_from="epoch_created",key_to="epoch_index")
        simulated_trees[chrom] = sort_tree_by_copy_number(simulated_trees[chrom])
        estimated_trees[chrom] = sort_tree_by_copy_number(estimated_trees[chrom])
        for i in range(5):
            pretty_print("##### CHROM: " + str(chrom))
        #pretty_print("The estimated tree is: " + str(estimated_trees[chrom]))
        #pretty_print("The simulated tree is: " + str(simulated_trees[chrom]))
        pretty_print("The simulated tree looks like: " + str(simulated_trees[chrom]))
        pretty_print("The simulated tree looks like: " + str(convert_dict_tree_to_list(simulated_trees[chrom])))
        pretty_print("The estimated tree looks like: " + str(estimated_trees[chrom]))
        pretty_print("The estimated tree looks like: " + str(convert_dict_tree_to_list(estimated_trees[chrom],total_epochs=pre+mid+post+(mid>=0)+(post>=0),is_truth=True)))
        if is_the_same_dict_tree_by_epoch_and_time_created(estimated_trees[chrom],simulated_trees[chrom]):
            pretty_print("They are the same")
        else:
            pretty_print("They are NOT the same")

        tree_sim = count_nodes_with_same_copy_number(estimated_trees[chrom], simulated_trees[chrom])    
        num_chrom_with_correct_CN += tree_sim #== sim_tree_length)
        sim_tree_len = count_nodes(simulated_trees[chrom])
        total_nodes += sim_tree_len
        pretty_print("In total there are " + str(tree_sim) + " nodes that have the exact same copy numbers out of " + 
                str(sim_tree_len))
        pretty_print("So this tree is " + str(int(float(tree_sim)*100/float(sim_tree_len))) + "% correct.")
        tree_sim = count_nodes_with_same_properties(estimated_trees[chrom], simulated_trees[chrom])    
        num_chrom_with_correct_CN_and_epoch_created += tree_sim #== tree_sim_length)

        pretty_print("In total there are " + str(tree_sim) + " nodes that have the exact same copy numbers and epochs created out of " + 
                str(sim_tree_len))
        pretty_print("So this tree is " + str(int(float(tree_sim)*100/float(sim_tree_len))) + "% correct.")

        pretty_print("It is quite hard to get the tree exactly correct so of the nodes that are correct by copynumber sum the difference in epoch_created")
        tree_sim = sum_tree_distance(estimated_trees[chrom], simulated_trees[chrom],diff_struct_is_inf = False)
        average_distance_from_truth_of_epoch_created += float(tree_sim)/float(sim_tree_len)/23
        pretty_print("The average distance per node is " + str(float(tree_sim)/float(sim_tree_len)))


    pretty_print("The total number of nodes in the tree is:")
    pretty_print(total_nodes)

    pretty_print("The number of nodes that match the true tree structure and have the correct CNs in our best estimate is:")
    pretty_print(num_chrom_with_correct_CN)

    pretty_print("The number of nodes that match the true tree structure and have the correct CNs AND the correct estimate for epoch created in our best estimate is:")
    pretty_print(num_chrom_with_correct_CN_and_epoch_created)

    pretty_print("The average distance of epoch created across all nodes to the true value is:")
    pretty_print(average_distance_from_truth_of_epoch_created)

    new_result = {"average_distance_from_truth_of_epoch_created":average_distance_from_truth_of_epoch_created,
            "num_chrom_with_correct_CN_and_epoch_created":num_chrom_with_correct_CN_and_epoch_created,
            "num_chrom_with_correct_CN":num_chrom_with_correct_CN,
            "total_nodes":total_nodes,
            "estimated_trees":estimated_trees,
            "simulated_trees":simulated_trees,
            "res":res,
            "val":val, 
            "pre_est":pre_est,
            "mid_est":mid_est, 
            "post_est":post_est, 
            "p_up_est":p_up_est, 
            "p_down_est":p_down_est, 
            "plambda_est":plambda_est, 
            "result":result,
            "BP_likelihoods":BP_likelihoods,
            "total_SNV_likelihood":total,
            "best_SNV_likelihood":best,
            "pre":real_pre,
            "mid":real_mid,
            "post":real_post,
            "p_up":real_p_up,
            "p_down":real_p_down,
            "plambda":real_rate
            }
    for key in new_result:
        pretty_print(key+": "+str(new_result[key]))

    d = shelve.open('file3_'+str(test_case)+'.txt')           
    if 'all_results' in d:
        all_results = d['all_results']

    all_results[str(result_index)] = new_result
    d['all_results'] = all_results
    d.close()




