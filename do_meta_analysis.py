import numpy as np
import glob
import os
import shelve
import sys


def convert_dict_tree_to_list(tree,total_epochs=None,is_truth=False):

    if tree is None:
        return None

    copy_number = tree.get('copy_number')
    epoch_created = tree.get('epoch_created')

    if is_truth:
        if 'child' in tree and tree["child"] is not None:
            epoch_created = tree["child"]["epoch_created"]
        else:
            epoch_created = total_epochs #max_epochs

    child_tree = convert_dict_tree_to_list(tree.get('child'),total_epochs,is_truth)
    complement_tree = convert_dict_tree_to_list(tree.get('complement'),total_epochs,is_truth)
    if child_tree is None:
        return [(copy_number, epoch_created)]
    return [(copy_number, epoch_created), child_tree, complement_tree]

def get_filenames_sorted_by_time():
    files = glob.glob("file3_*.txt.dir")

    wildcard_values = []

    for file in files:
        wildcard_value = os.path.splitext(os.path.basename(file))[0].split("_")[1].split(".")[0]
        creation_time = os.path.getctime(file)
        wildcard_values.append((wildcard_value, creation_time))

    wildcard_values = sorted(wildcard_values, key=lambda x: x[1])

    sorted_values = [val[0] for val in wildcard_values]

    return sorted_values

# Example usage:

def sort_dicts_by_aic(dicts):
    return sorted(dicts, key=lambda x: x['AIC'])

def sort_dicts_by_worst_aic(dicts):
    return sorted(dicts, key=lambda x: -x['AIC'])

def sort_dicts_by_val(dicts_list):
    return sorted(dicts_list, key=lambda x: x['val'])

def is_simulation_successful(simulation):
    mid_success = (simulation['mid_est'] != -1) == (simulation['mid'] != -1)
    post_success = (simulation['post_est'] != -1) == (simulation['post'] != -1)
    return mid_success and post_success

def get_successful_simulations(dicts_list):
    return [simulation for simulation in dicts_list if is_simulation_successful(simulation)]

def get_ev_string(pre, mid, post):
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

def count_genome_doublings(ev_string):
    return ev_string.count("G")

def add_ev_strings_and_counts_to_dicts(dicts_list):
    for d in dicts_list:
        d['ev_str'] = get_ev_string(d['pre'], d['mid'], d['post'])
        d['ev_str_est'] = get_ev_string(d['pre_est'], d['mid_est'], d['post_est'])
        d['genome_doublings'] = count_genome_doublings(d['ev_str'])
        d['genome_doublings_est'] = count_genome_doublings(d['ev_str_est'])
    return dicts_list

def create_3x3_grid(dicts_list):
    grid = np.zeros((3, 3))
    for d in dicts_list:
        if d == None:
            continue
        gd = d['genome_doublings']
        gd_est = d['genome_doublings_est']
        if gd < 3 and gd_est < 3:
            grid[gd, gd_est] += 1
    return grid

def find_matching_genome_doublings_index(dicts_list):
    for index, d in enumerate(dicts_list):
        if d['genome_doublings'] == d['genome_doublings_est']:
            return index
    return -1  # Indicates that no matching entry was found


def find_matching_ev_strings_index(dicts_list):
    for index, d in enumerate(dicts_list):
        if d['ev_str'] == d['ev_str_est']:
            return index
    return -1  # Indicates that no matching entry was found


def get_all_best_estimates(list_of_lists):
    return [x[0] for x in list_of_lists]
    
def get_all_true_estimates(list_of_lists):
    results = []
    for x in list_of_lists:
        index_matching_ev_string = find_matching_ev_strings_index(x)
        if index_matching_ev_string != -1:
            results += [list_of_lists[index_matching_ev_string]]
        else:
            results += [None]
    return results

def add_aic_to_dicts(dicts):
    for d in dicts:
        num_parameters = d['genome_doublings_est'] + 3
        neg_log_likelihood = d['val']
        aic = 2 * num_parameters + 2 * neg_log_likelihood
        d['AIC'] = aic
    return dicts 


import shelve

def process_test_cases(sorted_values):
    count = 0
    list_of_lists = []

    for test_case in sorted_values:
        print(test_case)

        d = shelve.open('file3_' + str(test_case) + '.txt')
        if len(list(d.keys())) > 0:
            all_results = d['all_results']

            print("length of results:" + str(len(all_results.keys())))
            all_results = [all_results[x] for x in all_results]
            for result in all_results:
                result['test_case'] = test_case

            all_results = add_ev_strings_and_counts_to_dicts(all_results)
            all_results = add_aic_to_dicts(all_results)
            all_results = sort_dicts_by_aic(all_results)

            for result in all_results:
                print_result_info(result)

            count += 1

        list_of_lists += [all_results]
        d.close()

    return list_of_lists

def print_result_info(result):
    fields = ["pre", "mid", "post", "p_up", "p_down"]
    est_fields = [f + "_est" for f in fields]

    print("(truth,estimated)")
    for field, est_field in zip(fields, est_fields):
        print(f"{field}: {result[field], result[est_field]}")

    print("val:" + str(result["val"]))
    print("AIC:" + str(result["AIC"]))
    print("ev_string:" + str(result["ev_str"]))


sorted_filenames = get_filenames_sorted_by_time()
list_of_lists = process_test_cases(sorted_filenames)

def top_n_lists_by_best_aic(lists_of_dicts, n):
    sorted_lists = sorted(lists_of_dicts, key=lambda l: min(d['AIC'] for d in l))
    return sorted_lists[:n]


import random

def filter_lists_by_ev_str(lists_of_dicts, p):
    filtered_lists = [
            list_of_dicts for list_of_dicts in lists_of_dicts
            if list_of_dicts[0]['ev_str_est'] == list_of_dicts[0]['ev_str'] or random.random() < p 
    ]
    return filtered_lists
    #return [l for l in filtered_lists if l]  # Remove empty lists

# Example usage
#filtered_lists = filter_lists_by_ev_str(your_lists_of_dicts, 0.5)

print("RAW")
all_true_estimates = get_all_true_estimates(list_of_lists)
all_best_estimates = get_all_best_estimates(list_of_lists)

worst_are_first = sort_dicts_by_worst_aic(all_best_estimates)



sys.exit()


print("BEST ESTIMATED")
print(create_3x3_grid(all_best_estimates))
print("IF ONLY TRUTH EXISTED")
print(create_3x3_grid(all_true_estimates))


list_of_lists = filter_lists_by_ev_str(list_of_lists,0.2)
all_true_estimates = get_all_true_estimates(list_of_lists)
all_best_estimates = get_all_best_estimates(list_of_lists)

print("BEST ESTIMATED")
print(create_3x3_grid(all_best_estimates))
print("IF ONLY TRUTH EXISTED")
print(create_3x3_grid(all_true_estimates))

GD_indexes = [find_matching_genome_doublings_index(x) for x in list_of_lists]
path_indexes = [find_matching_ev_strings_index(x) for x in list_of_lists]

