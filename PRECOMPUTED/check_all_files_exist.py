import os

all_files_exist = True

for x in range(101):
    for y in range(101-x):
        subbed_mat_base = f"subbed_mat_u{x}_d{y}_p8_v4"
        collated_base = f"collated_u{x}_d{y}_p8_v4"
        subbed_mat_types = [".pickle", ".powers.pickle", ".precomputed_paths.pickle"]
        collated_types = [".csv"]
        missing = []
        for ext in subbed_mat_types:
            filename = os.path.join("MATRICES", subbed_mat_base + ext)
            if not os.path.isfile(filename):
                missing.append(ext)
        for ext in collated_types:
            filename = os.path.join("MATRICES", collated_base + ext)
            if not os.path.isfile(filename):
                missing.append(ext)
        if missing:
            all_files_exist = False
            print(f"Missing files for x={x}, y={y}: {', '.join(missing)}")

if all_files_exist:
    print("All files exist for every combination of x and y.")

