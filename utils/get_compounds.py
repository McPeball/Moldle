import pubchempy as pcp

with open("../resources/compounds10k.tsv", "w") as f:
    #f.write(f'iupac_name\tisomeric_smiles\n')
    for i in range(1, 10000):
        print(i)
        c = pcp.Compound.from_cid(i)
        f.write(f'{c.iupac_name}\t{c.isomeric_smiles}\n')