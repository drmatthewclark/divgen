from divgen import makexample, sdfile
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import random

EXPAND_ISOMERS = False
target = int(1e2)
count = 0
mset = set() 

while True:
    try:
        mol = makexample(random.randrange(10,50))
        sd = sdfile(mol)
        molecule  = Chem.MolFromMolBlock(sd, sanitize = True)
        firstsmiles = Chem.MolToSmiles(molecule, canonical = True, isomericSmiles = True)
        shash = hash(firstsmiles)

        if shash in mset:
            continue 

        if EXPAND_ISOMERS:
            print(firstsmiles) 
            isomers = EnumerateStereoisomers(molecule)
        else:
            isomers = [molecule]
    
        for rdkit_mol in isomers:
            smiles = Chem.MolToSmiles(rdkit_mol, canonical = True, isomericSmiles = True)
            mset.add(hash(smiles))
            print(smiles)

            count += 1
            if count >=  target:
                exit(0)         

    except Exception as e:
        #print(e)
        pass

