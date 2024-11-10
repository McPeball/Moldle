import editdistance
import random

def read_compounds(file = "../resources/compounds.tsv"):
    with open(file, "r") as f:
        compound_names = []
        smiles = []
        for l in f:
            data_elements = l.split("\t")
            if len(data_elements[0]) > 30:
               continue
            compound_names.append(data_elements[0])
            smiles.append(data_elements[1])
    compounds = {
        "name":compound_names,
        "SMILES":smiles
    }
    return compounds

# pick a compound at random from the dataset
def choose_compound_index(compounds):
    compound_index = random.randint(0, (len(compounds["name"]) - 1))
    
    return compound_index

# Compute the edit distance between the random molecule and all others
def compute_edit_dists(compounds, compound):
    dists = dict()
    for c in compounds['name']:
        dists[c] = editdistance.eval(c, compound)
        #dists.append(editdistance.eval(c, compound))
    return dists

def main():
    compounds = read_compounds()
    compound_index = choose_compound_index(compounds)
    
    # compute Levenstein (edit) distances between compound names
    dists = compute_edit_dists(compounds, compounds['name'][compound_index])
    
    # sort the dists dict by values
    dists = dict(sorted(dists.items(), key=lambda item: item[1]))
    
    # print the five closest compound names
    names = list(dists.keys())
    for i in range(0, 5):
        print(names[i])

if __name__ == '__main__':
    main()




