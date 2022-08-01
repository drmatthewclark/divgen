from divgen import makexample, sdfile
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import random
target = int(1e2)
count = 0
mset = set() 

while True:
    try:
        mol = makexample(random.randrange(10,50))
        sd = sdfile(mol)
        rdkit_mol  = Chem.MolFromMolBlock(sd, sanitize = True)
        firstsmiles = Chem.MolToSmiles(rdkit_mol, canonical = True, isomericSmiles = True)
        shash = hash(firstsmiles)

        if shash in mset:
            continue 

        print(firstsmiles) # no isomer version

        #isomers = EnumerateStereoisomers(rdkit_mol)
        isomers = rdkit_mol
        for rdkit_mol in isomers:
            smiles = Chem.MolToSmiles(rdkit_mol, canonical = True, isomericSmiles = True)
            mset.add(hash(smiles))
            print(smiles)

            count += 1
            if count ==  target:
                    break

    except Exception as e:
        #print(e)
        pass

