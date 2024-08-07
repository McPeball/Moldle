import streamlit as st
import random
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from scipy import spatial as sc
#import scipy.spatial.distance.euclidean as euc

# Get dataset of smiles and iupac names
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

# show the compound using rdkit

def draw_smiles(compounds, compound_index, png_file = "compound.png") -> None:
    m = Chem.MolFromSmiles(compounds["SMILES"][compound_index])
    img = Draw.MolToImage(m)
    img.save(png_file)
    st.image(png_file)

# find plausible similar names (red herrings)
# break name into parts, compare other compound name parts
# score based on count of common parts
def get_red_herrings(compounds, compound_index):
    red_herrings = set()
    while len(red_herrings) < 4:
        j = random.randint(0, (len(compounds["SMILES"]) - 1))
        if j == compound_index:
            continue
        red_herrings.add(j)
    red_herrings = list(red_herrings)
    red_herrings.append(compound_index)
    red_herrings.sort()
    return red_herrings

# Another red-herring idea:
# Read the full set of iupac_names and break each into parts
# Log the relationship between compound names based on their component
# eg. If 'pyro' is component of the target compound name, it and any other
# compounds that also contain 'pyro' would get a count of one for 'pyro'
# Other compounds would get a count of zero.








# print correct name and some red herrings, ask user
# to guess the corrent name
def user_guess(compounds, red_herrings):
    compound_guess = st.radio(
        "**Select a compound**",
        [
            compounds["name"][red_herrings[0]],
            compounds["name"][red_herrings[1]],
            compounds["name"][red_herrings[2]],
            compounds["name"][red_herrings[3]],
            compounds["name"][red_herrings[4]]
        ],
        index=None,
    )

def check_user_guess(compound_guess, compounds, compound_index):
    #st.write(compound_guess)
    #st.write(compounds["name"[compound_index]])
    if str(compound_guess) == compounds["name"][compound_index]:
        st.write("yes!")
    else:
        st.write("no!")

def main():
    compounds = read_compounds()
    print(compounds)
    compound_index = choose_compound_index(compounds)
    draw_smiles(compounds, compound_index, png_file = "compound.png")
    red_herrings = get_red_herrings(compounds, compound_index)
    compound_guess = None
    compound_guess = user_guess(compounds, red_herrings)
    check_user_guess(compound_guess, compounds, compound_index)

if __name__ == "__main__":
    main()