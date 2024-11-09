import streamlit as st
import random
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from scipy import spatial as sc

# Get dataset of smiles and iupac names
def read_compounds(file = "../resources/compounds.tsv", max_name_length = 100):
    """Read a tsv file of compound names and isomeric SMILES strings
    The input file is generated using utils/get_compounds.py, reading
    data from Pubchem.

    Returns a dict with two keys: 'name' and 'SMILES'. Each key is associated
    with a list of IUPAC names or isomeric SMILES strings. The two lists of values
    correpsond with each other. ie. The first element of the list 'name' corresponds
    to the first element of the list in 'SMILES' and so on.
    """
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


def choose_compound(compounds) -> None:
    """Picks a random index within the range of available compounds. It then uses this
    index to select a name and corresponding SMILES string and store these in the
    st.session_state dict with keys 'compound_name' and 'compound_smiles'.
    """
    compound_index = random.randint(0, (len(compounds["name"]) - 1))
    st.session_state.compound_name = st.session_state.compounds['name'][compound_index]
    st.session_state.compound_smiles = st.session_state.compounds['SMILES'][compound_index]


def draw_smiles(compounds, compound_name, png_file = "compound.png") -> None:
    """Draws a chemical structure using the current string stored in
    st.session_state.compound_smiles and saves the image to a file called
    compound.png. The png file is then rendered to the page.
    """
    m = Chem.MolFromSmiles(st.session_state.compound_smiles)
    img = Draw.MolToImage(m, size=(300, 150))
    img.save(png_file)
    st.image(png_file)

def get_red_herrings(compounds, compound_name):
    """Selects a random set of five compound names, including the target compound name.
    Returns a list of the selected compound names sorted alphabetically, to randomise the
    position of the target compound.
    """
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
    """Checks if the user guess matches the target compound name.
    If the match is correct, the score gets incremented. Otherwise it is reset to zero
    and the 'wrong' and 'should_have_been' values are set to be displayed later.
    Finally, a new compound name and corresponding SMILES string is selected.
    """
    if st.session_state.guess == st.session_state.compound_name:
        st.session_state.user_score += 1
    else:
        st.session_state.wrong = True
        st.session_state.should_have_been = st.session_state.compound_name
        st.session_state.user_score = 0
    choose_compound(st.session_state.compounds)
    st.session_state.red_herrings = get_red_herrings(st.session_state.compounds, st.session_state.compound_name)
    

def user_guess(compounds, red_herrings, compound_index):
    """Creates a form with a radio button set to select a compound name
    and a submit button to call the check_user_guess() function.
    """
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
    """Main function to create the page.
    """
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