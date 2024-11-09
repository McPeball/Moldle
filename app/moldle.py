import streamlit as st
import random
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from scipy import spatial as sc

# Get dataset of smiles and iupac names
def read_compounds(file = "../resources/compounds.tsv", max_name_length = 100):
    with open(file, "r") as f:
        compound_names = []
        smiles = []
        for l in f:
            data_elements = l.split("\t")
            if len(data_elements[0]) > max_name_length:
                continue
            compound_names.append(data_elements[0])
            smiles.append(data_elements[1])
    compounds = {
        "name":compound_names,
        "SMILES":smiles
    }
    return compounds

# pick a compound at random from the dataset
def choose_compound(compounds) -> None:
    compound_index = random.randint(0, (len(compounds["name"]) - 1))
    st.session_state.compound_name = st.session_state.compounds['name'][compound_index]
    st.session_state.compound_smiles = st.session_state.compounds['SMILES'][compound_index]

# show the compound using rdkit
def draw_smiles(compounds, compound_name, png_file = "compound.png") -> None:
    m = Chem.MolFromSmiles(st.session_state.compound_smiles)
    img = Draw.MolToImage(m, size=(300, 150))
    img.save(png_file)
    st.image(png_file)

# find plausible similar names (red herrings)
# break name into parts, compare other compound name parts
# score based on count of common parts
def get_red_herrings(compounds, compound_name):
    red_herrings = set([st.session_state.compound_name])
    while len(red_herrings) < 5:
        rh_index = random.randint(0, (len(compounds['name']) - 1))
        rh_name = st.session_state.compounds['name'][rh_index]
        if rh_name == st.session_state.compound_name:
            continue
        red_herrings.add(rh_name)
    red_herrings = list(red_herrings)
    red_herrings.sort()
    return red_herrings

def check_user_guess() -> None:
    if st.session_state.guess == st.session_state.compound_name:
        st.session_state.user_score += 1
    else:
        st.session_state.wrong = True
        st.session_state.should_have_been = st.session_state.compound_name
        st.session_state.user_score = 0
    choose_compound(st.session_state.compounds)
    st.session_state.red_herrings = get_red_herrings(st.session_state.compounds, st.session_state.compound_name)
    

def user_guess(compounds, red_herrings, compound_index):
    st.session_state.guess = st.radio(
        "**Choose a name:**",
        st.session_state.red_herrings,
        index=None
    )
    check_guess = st.button(
        "Submit",
        on_click = check_user_guess
        )

def main():
    if "wrong" not in st.session_state:
        st.session_state.wrong = False
        st.session_state.should_have_been = ''
    if "user_score" not in st.session_state:
        st.session_state.user_score = 0
    if "compounds" not in st.session_state:
        st.session_state.compounds = read_compounds()
    if 'compound_name' not in st.session_state:
        choose_compound(st.session_state.compounds)
    if 'red_herrings' not in st.session_state:
        st.session_state.red_herrings = get_red_herrings(st.session_state.compounds, st.session_state.compound_name)
    st.image('../resources/logo.png', width=300)
    st.write('**Name this molecule:**')
    draw_smiles(st.session_state.compounds, st.session_state.compound_name, png_file = "compound.png")
    user_guess(st.session_state.compounds, st.session_state.red_herrings, st.session_state.compound_name)
    if st.session_state.wrong:
        st.write(f"The correct answer was {st.session_state.should_have_been}")
        st.session_state.wrong = False
    st.write(f"**Streak:** {st.session_state.user_score}")

if __name__ == "__main__":
    main()