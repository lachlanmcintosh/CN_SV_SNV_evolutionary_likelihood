
import glob
import os
import shelve

files = glob.glob("../file3_*.txt.dir")

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
    
    d = shelve.open('../file3_'+str(test_case)+'.txt')
    if len(list(d.keys())) > 0:
        all_results = d['all_results']
        print("length of results:"+str(len(all_results.keys())))
        for key in all_results:
            print("key:" + str(key))
            print(all_results[key])
        count += 1

        if count >= 10:
            break

    d.close()


