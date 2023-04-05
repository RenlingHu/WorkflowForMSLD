import os

def write_mol_list(filename):
    with open ("./mol_list.txt","a+") as fp:
        fp.write("{}\n".format(filename))

# clean previous content in mol_list.txt
with open ("./mol_list.txt","a+") as fp0:
    fp0.seek(0)
    fp0.truncate()

# read in the name of mol2_file in work directory
if __name__ == "__main__":
    fp = os.walk(".")
    for path, dir_list, file_list in fp:
        for file_name in file_list:
            if ".mol2" in file_name:
                write_mol_list(file_name[:-5])