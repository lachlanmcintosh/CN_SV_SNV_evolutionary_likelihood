
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


files = glob.glob("file3_*.txt.dir")

wildcard_values = []

for file in files:
    wildcard_value = os.path.splitext(os.path.basename(file))[0].split("_")[1].split(".")[0]
    creation_time = os.path.getctime(file)
    wildcard_values.append((wildcard_value, creation_time))

wildcard_values = sorted(wildcard_values, key=lambda x: x[1])

sorted_values = [val[0] for val in wildcard_values]

print(sorted_values)

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


def get_all_best_estiamtes(list_of_lists):
    return [x[0] for x in list_of_lists]
    
def get_all_true_estiamtes(list_of_lists):
    results = []
    for x in list of lists:
        index_matching_ev_string = find_matching_ev_strings_index(x)
        if index_matching_ev_string != -1:
            results += x[index_matching_ev_string]
        else:
            results += None
 
    return results

def add_aic_to_dicts(dicts, num_parameters):
    for d in dicts:
        num_parameters = d['genome_doublings_est'] + 3
        neg_log_likelihood = d['val']
        aic = 2 * num_parameters + 2 * neg_log_likelihood
        d['AIC'] = aic
    return dicts 


count = 0
for test_case in sorted_values:
    print(test_case)
    
    d = shelve.open('file3_'+str(test_case)+'.txt')
    print(d.keys())
    if len(list(d.keys())) > 0:
        all_results = d['all_results']
        print("length of results:"+str(len(all_results.keys())))
        
        #key = "0"
        #if True: #all_results[key]["val"] != "inf": 
        for key in all_results:
            print("key:" + str(key))
            #print(all_results[key])

            pre = all_results[key]["pre"]
            mid = all_results[key]["mid"]
            post = all_results[key]["post"]
            p_up = all_results[key]["p_up"]
            p_down = all_results[key]["p_down"]
            #plambda = all_results[key]["plambda"]

            pre_est = all_results[key]["pre_est"]
            mid_est = all_results[key]["mid_est"]
            post_est = all_results[key]["post_est"]
            p_up_est = all_results[key]["p_up_est"]
            p_down_est = all_results[key]["p_down_est"]
            #plambda_est = all_results[key]["plambda_est"]

            print("(truth,estimated)")
            print("pre:"+str((pre,pre_est)))
            print("mid:"+str((mid,mid_est)))
            print("post:"+str((post,post_est)))
            print("p_up:"+str((p_up,p_up_est)))
            print("p_down:"+str((p_down,p_down_est)))
            #print("plambda:"+str((plambda,plambda_est)))

            val = all_results[key]["val"]
            print("val:"+str(val))


            total_epochs = all_results[key]["total_epochs"] #pre+mid+post+(pre>=0)+(post>=0)


            est = all_results[key]["estimated_trees"]
            sim = all_results[key]["simulated_trees"]
            #for chrom in est:
            #    print(convert_dict_tree_to_list(est[chrom]))
            #    print(convert_dict_tree_to_list(sim[chrom]))

        count += 1

        if count >= 10:
            break

    d.close()

#for test_case in sorted_values:
#    print(test_case)
#    d = shelve.open('file3_'+str(test_case)+'.txt')
#    if len(list(d.keys())) > 0:
#        all_results = d['all_results']
#        val_key_pairs = [(all_results[key]["val"],key) for key in all_results] 
#
#        best_key = sorted(val_key_pairs)[0][1]
	


	

