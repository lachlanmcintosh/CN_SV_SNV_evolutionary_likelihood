
import glob
import os
import shelve


def convert_dict_tree_to_list(tree):
    if tree is None:
        return None

    copy_number = tree.get('copy_number')
    epoch_created = tree.get('epoch_created')

    if 'child' in tree and tree["child"] is not None:
        epoch_killed = tree["child"]["epoch_created"]
    else:
        epoch_killed = total_epochs #max_epochs

    child_tree = convert_dict_tree_to_list(tree.get('child'))
    complement_tree = convert_dict_tree_to_list(tree.get('complement'))
    if child_tree is None:
        return [(copy_number, epoch_killed)]
    #return [(copy_number, epoch_created), child_tree, complement_tree]
    return [(copy_number, epoch_killed), child_tree, complement_tree]



files = glob.glob("file3_*.txt.dir")

wildcard_values = []

for file in files:
    wildcard_value = os.path.splitext(os.path.basename(file))[0].split("_")[1].split(".")[0]
    creation_time = os.path.getctime(file)
    wildcard_values.append((wildcard_value, creation_time))

wildcard_values = sorted(wildcard_values, key=lambda x: x[1])

sorted_values = [val[0] for val in wildcard_values]

print(sorted_values)

count = 0
for test_case in sorted_values:
    print(test_case)
    
    d = shelve.open('file3_'+str(test_case)+'.txt')
    if len(list(d.keys())) > 0:
        all_results = d['all_results']
        print("length of results:"+str(len(all_results.keys())))
        for key in all_results:
            print("key:" + str(key))
            print(all_results[key])
            est = convert_dict_tree_to_list(all_results[key]["estimated_trees"])
            sim = convert_dict_tree_to_list(all_results[key]["simulated_trees"])
            for chrom in est:
                print("chrom:"+str(chrom))
                print("simulated tree:"+str(sim[chrom]))
                print("estimated tree:"+str(est[chrom]))

        count += 1

        if count >= 10:
            break

    d.close()


