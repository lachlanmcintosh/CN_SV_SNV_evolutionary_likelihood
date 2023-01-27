
length = 128
all_paths = []
for i in range(129):
    all_paths.append(str(i))

G1 = list(range(0,129))# + [128]

for i in range(len(G1)):
    for j in range(len(G1)):
        if G1[i] + G1[j] <= length:
            all_paths.append(str(G1[i]) + "G" + str(G1[j])) 

G2 = list(range(0,9))+list(range(10,17,2))+list(range(20,33,4)) + list(range(40,65,8)) + list(range(80,129,16)) #[32,64,128] #list(range(6,12,2))+list(range(20,110,10))
       
print(G2)
 
for i in range(len(G2)):
    for j in range(len(G2)):
        for k in range(len(G2)):
            if G2[i] + G2[j] + G2[k] <= length:
                all_paths.append(str(G2[i]) + "G" + str(G2[j]) + "G" + str(G2[k])) 

G2 = set(G2)

print(G2)
print(len([x for x in all_paths if len(x.split("G")) == 1]))
print(len([x for x in all_paths if len(x.split("G")) == 2]))
print(len([x for x in all_paths if len(x.split("G")) == 3]))
print(len([x for x in all_paths]))

# save the like of all combinations that we wantto make

save(all_paths,"all_path_combinations_p128_v3")
save([G1,G2],"all_path_combinations_p128_v3_lists")

