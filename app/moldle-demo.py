import streamlit as st
import random
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from scipy import spatial as sc
#import scipy.spatial.distance.euclidean as euc

# Get dataset of smiles and iupac names
with open("compounds.tsv", "r") as f:
    compounds = []
    smiles = []
    for l in f:
        data_elements = l.split("\t")
        if len(data_elements[0]) > 30:
            continue
        compounds.append(data_elements[0])
        smiles.append(data_elements[1])

# pick a compound at random from the dataset
i = random.randint(0, (len(compounds) - 1))

# show the compound using rdkit
print(smiles[i])

m = Chem.MolFromSmiles(smiles[i])
img = Draw.MolToImage(m)
img.save("compound.png")
st.image("compound.png")

# find plausible similar names (red herrings)
# break name into parts, compare other compound name parts
# score based on count of common parts

red_herrings = set()
while len(red_herrings) < 4:
    j = random.randint(0, (len(compounds) - 1))
    if j == i:
        continue
    red_herrings.add(j)
red_herrings = list(red_herrings)
red_herrings.append(i)
red_herrings.sort()
print(red_herrings)
print(i)

# Another red-herring idea:
# Read the full set of iupac_names and break each into parts
# Log the relationship between compound names based on their component
# eg. If 'pyro' is component of the target compound name, it and any other
# compounds that also contain 'pyro' would get a count of one for 'pyro'
# Other compounds would get a count of zero.








# print correct name and some red herrings, ask user
# to guess the corrent name

compound_guess = st.radio(
    "**Select a compound**",
    [compounds[red_herrings[0]], compounds[red_herrings[1]], compounds[red_herrings[2]], compounds[red_herrings[3]], compounds[red_herrings[4]]],
    index=None,
)

st.write(compound_guess)
st.write(compounds[i])
if str(compound_guess) == compounds[i]:
    st.write("yes!")
else:
    st.write("no!")