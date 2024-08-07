import streamlit as st

some_names = ["Alice", "Bob2", "Charlie"]

genre = st.radio(
    "What's your favorite movie genre",
    some_names,
)

if genre == "Alice":
    st.write("You selected Alice.")
else:
    st.write("You didn't select Alice.")