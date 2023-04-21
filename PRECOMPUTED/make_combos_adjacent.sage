import sys

length = int(sys.argv[1]) 
base_paths = list(range(0,length+1))
all_paths = []
for i in base_paths:
    all_paths.append(str(i))

for i in base_paths:
    for j in base_paths:
        if i + j + 1 <= length:
            all_paths.append(str(i) + "G" + str(j)) 

for i in base_paths:
    for j in base_paths:
        for k in base_paths:
            if i + j + k + 2 <= length:
                all_paths.append(str(i) + "G" + str(j) + "G" + str(k)) 

print(base_paths)
print(len([x for x in all_paths if len(x.split("G")) == 1]))
print(len([x for x in all_paths if len(x.split("G")) == 2]))
print(len([x for x in all_paths if len(x.split("G")) == 3]))
print(len([x for x in all_paths]))

# save the like of all combinations that we wantto make

save(all_paths,"all_path_combinations_p"+str(length)+"_v4")

