#!/usr/bin/env python
# coding: utf-8

# Divgen  Matthew Clark  2022
# CC BY-NC-SA 4.0
# https://creativecommons.org/licenses/by-nc-sa/4.0/

# In[1]:


import random
import os
import math
from datetime import datetime
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit import Chem
import time

# bond - tuple (fromatom, toatom, bondorder, chiral)  sorted so that toatom > fromatom
# atom - ordered list of atom names, e.g. 'C'
# molecule - tuple of  (atom, bond)
ringatoms = ['C', 'O', 'S', 'N', 'P' ]
atnum_map = {'H' : 1, 'Li' : 3, 'Be' : 4, 'B' : 5, 'C' : 6, 'N' : 7, 'O' : 8, 'F' : 9, 'Mg' : 12, 'Si' : 14, 'Cl' : 32, 'S' : 16, 'P' : 15 , 'Br' : 35, 'I' : 53}
add_hetero = True
PROB_AROMATIC = 0.3
PROB_NAPHTHYL =  0.1
PROB_SECOND_AROMATIC = 0.1
PROB_CARBOXYLIC = 0.1
GEN_ISOMERS = True


# hetero atoms, valence, and probabilities for them
hetatoms = {}
hetatoms['S']  =  (2,0.200)
hetatoms['Br'] =  (1,0.05)
hetatoms['I']  =  (1,0.05)
hetatoms['Cl'] =  (1,0.154)
hetatoms['F']  =  (1,0.167)
hetatoms['N']  =  (3,0.645)
hetatoms['O']  =  (2,0.600)
hetatoms['P']  =  (3,0.050)
hetatoms['C']  =  (4,0.850)

def counthet(atoms):
    count = 0
    for s in atoms:
        if s != 'C': count += 1
        
    return count

def chooseHetero(nhet):
    """
    choose a hetero replacement based on probabilities in the
    map
    """
    a = max(1.0,nhet*nhet)
 
    elements = ['S'    , 'Br',   'I',   'Cl',    'F',   'N',     'O',   'P' , 'C'  ]
    weights =  [0.137/a, .02/a, .05/a,  0.05/a, 0.05/a, .545/a, .644/a, .050/a , 0.98]
    
    symbol =  random.choices(elements, weights = weights, k = 1)[0]
    (valence, prob) = hetatoms[symbol]

    return (symbol, hetatoms[symbol])

def heteroatoms(molecule):
    """
    randomly assign some carbons to be heteroatoms based on
    map of probabilities
    """
    atoms, bonds = molecule
     
    # walk over atoms in random order, choose an hetero
    for atom in random.sample(range(len(atoms)), len(atoms)):
        
        if atoms[atom] != 'C': # don't re-change
            continue
            
        het = chooseHetero(counthet(atoms))
        
        symbol, (val, prob) = het
        if symbol == 'C':
            continue
        current_valence = valence(bonds, atom)

        if current_valence <=  val:
            atoms[atom] = symbol

    return (atoms, bonds)


# In[3]:


def valence(bonds, atnum):
    """
        count the number of bonds to the given atom
    """
    valence = 0
    for (f, t, bondorder, chiral) in bonds:
        if f == atnum or t == atnum:
            valence += bondorder
    
    return valence


# In[4]:


def sdfile(molecule):
    """
    create a string from the molecule that is the SDFile (cf. CTFILE Formats) for this molecule
    """
    atoms, bonds = molecule
    bonds.sort()
    
    # one could assign a name and comment
    name = ''
    comment = ''
    
    result = "%-78s\n" % (name)
    userstring = '  MCMolgen'
    datestring = datetime.now().strftime('%m%d%y%H%M')
    result += userstring + datestring + '2D' + '\n'
    result += '%-80s\n' % (comment)
    result += "%3d%3d  0  0  0  0            999 V2000\n" % (len(atoms), len(bonds))
    
    for symbol in atoms:
        result += "%10.4f%10.4f%10.4f %3s               \n" % ( 0.0, 0.0, 0.0, symbol)
    
    for (f, t, bo, c) in bonds:
        result += "%3d%3d%3d%3d\n" % (f + 1, t + 1, bo, c)
    
    result += "M  END\n\n$$$$"
    
    return result


# In[12]:


def multiplebonds(molecule, nbond, border):
    """
    randomly select eligible bonds to upgrade to the bond order given by "border"
    """
    
    if nbond == 0:
        return molecule
   
    atoms, bonds = molecule
    
    nbond = min(len(bonds), nbond)
    
    adjust = 0
    if border == 3: adjust = -1
    
    # select nbond bonds randomly from all the bonds
    for index in random.sample(range(len(bonds)), nbond):
        
        f,t,bo,c = bonds[index]
         
        if bo == border: # already has this bondorder
            continue
            
        fv = valence(bonds, f)  # from atom valence
        tv = valence(bonds, t)  # to   atom valence
        
        if  4 - fv  >=  border + adjust         and 4 - tv  >=  border + adjust          and bo != border:
            # upgrade bond
            bonds[index] = (f, t, border, c)

    return (atoms, bonds)



def testbond(bonds, bond):
    """ 
        test that the bond is valid, that it does not
        already exist in any bond order, and that the atoms
        are not equal
    """
    
    bf,bt,o,c = bond
    
    if (bf == bt):
        return False
    
    for (f,t,o,c) in bonds:
        if (bf, bt) == (f,t):
            return False
        
    return True

    
def addbond(bonds, bond):
    """
    add the bond to the list of bonds
    """
    if testbond(bonds, bond):    
        bonds.append(bond)
        
    return bonds


def fixbond( bond ):
    """
    sort bond order so that the fromatom > toatom
    """
    a, b, c, d = bond
    
    if a > b:
        return (b, a, c, d)
    
    return (a, b, c, d)


    
def tracepath(molecule, start, path):
    """
    follow graph to create a connected path
    """
    atoms, bonds = molecule
    
    for (f, t, o, c) in random.sample(bonds, len(bonds)):
        if t == start:
            if not f in path:
                path.append(f)
                return tracepath(molecule, f, path)
                   
        if f == start:
            if not t in path:
                path.append(t)
                return tracepath(molecule, t, path)
            
    return (molecule, start, path)
 

    
def paths(molecule, start):
    """
    find some paths to create rings from the current atom
    this allows biasing for ring size
    """
    result = []
    s = start
    for i in range(3):
        path = [s]
        mol, start, test = tracepath(molecule, start, path)
        if len(test) > len(result):
            result = test
    
    return result
    
            
    
def ringbonds(molecule, nrings):
    """
    add nrings random ring closures
    """
    
    if nrings == 0:
        return molecule
    
    atoms, bonds = molecule
    
    for i in range(nrings):
        bondfrom = random.sample(range(len(atoms)), 1)[0]
        path = paths(molecule, bondfrom)
        if len(path) < 4:
            continue
       
        # pick a path to make a ring of average size 6
        bondto = path[min(gaussrandom(6, 1), len(path)-1)]
       
        if atoms[bondto] not in ringatoms or atoms[bondfrom] not in ringatoms:
            continue

        b = fixbond((bondfrom, bondto, 1, 0))

        if testbond(bonds, b)                     and valence(bonds, bondfrom) < 4          and valence(bonds, bondto)   < 4:
            bonds = addbond(bonds, b)

    return(atoms, bonds)



def makearomatic(molecule):
    """
    create a benzene ring
    """
    
    atoms, bonds = molecule
    
    bondto = random.randrange(0,len(atoms)) # pick an atom to attach to
    if valence(bonds, bondto) > 3 or  atoms[bondto] not in ringatoms:
        return molecule
    
    startatom = len(atoms)
    atoms.append('C')
    bond = fixbond( (startatom, bondto, 1, 0) ) 
    if not testbond(bonds, bond):
        print("problem with makearomatic bond")
        print(bond)
        
    bonds = addbond(bonds, bond)

    for i in range(startatom, startatom + 5):
        atoms.append('C')
        bond = (i, i + 1, (i%2) + 1, 0)
        bonds = addbond(bonds, bond)
    
    i += 1
    bond = (startatom, startatom + 5, (i%2) + 1, 0)
    bonds = addbond(bonds, bond)
   
    return (atoms, bonds)


def makenaphthyl(molecule):
    """
    create a naphthyl ring because it is statistically less
    likely than the frequency in molecules
    """
    
    atoms, bonds = molecule
    bondto = random.randrange(0,len(atoms)) # pick an atom to attach to
    if valence(bonds, bondto) > 3 or  atoms[bondto] not in ringatoms:
        return molecule

    startatom = len(atoms)
    atoms.append('C')
    bond = fixbond( (startatom, bondto, 1, 0) ) 
    if not testbond(bonds, bond):
        print("problem with makearomatic bond")
        print(bond)

    bonds = addbond(bonds, bond)

    for i in range(startatom, startatom + 5):
        atoms.append('C')
        bond = fixbond((i, i + 1, (i%2) + 1, 0))
        bonds = addbond(bonds, bond)
    
    i += 1
    bond = fixbond((startatom, startatom + 5, (i%2) + 1, 0))
    bonds = addbond(bonds, bond)
 
    for i in range(4):
        atoms.append('C')
        bond = fixbond((len(atoms)-2, len(atoms)-1, (i%2)+1, 0 ))
        bonds = addbond(bonds, bond)

    bond =  fixbond((len(atoms)-6, len(atoms)-1, 1, 0 ))
    bonds = addbond(bonds, bond)

    return (atoms, bonds)



def makecarbonyl(molecule):
    """
    create a carbonyl group.
    """
    atoms, bonds = molecule
    probability =  0.1
    
    # walk over atoms in random order, choose an hetero
    for atom in random.sample(range(len(atoms)), int(len(atoms)/2) + 1 ):
        rand = random.random() 
        if rand > probability:
            continue

        if atoms[atom] != 'C': # don't alter
            continue
            
        current_valence = valence(bonds, atom)
        
        # reduce probability of aldehydes
        if current_valence == 1 and rand  > probability:
            continue
            
        if current_valence < 3:
            atoms.append('O')
            bond =  fixbond((atom, len(atoms)-1, 2, 0))
            bonds = addbond(bonds, bond)
            
    return (atoms, bonds)


def makecarboxylic(molecule):
    """
    create a carboxylic acid group
    """
    
    atoms, bonds = molecule
    
    bondto = random.randrange(0,len(atoms)) # pick an atom to attach to
    
    if valence(bonds, bondto) > 3 or atoms[bondto] != 'C':
        return molecule
    
    newatom = len(atoms) # next atom #
    atoms.append('C')
    bond= fixbond(( bondto, newatom, 1, 0))
    bonds = addbond(bonds, bond)
                  
    newatom2 = len(atoms) # next atom #
    atoms.append('O')
    bond= fixbond(( newatom2, newatom, 2, 0))
    bonds = addbond(bonds, bond)                  
                  
    newatom3 = len(atoms) # next atom #
    atoms.append('O')
    bond= fixbond(( newatom, newatom3, 1, 0))
    bonds = addbond(bonds, bond)                            
                  
    return (atoms, bonds)


# In[11]:


def smilesprinter(rdkit_mol, fname = ''):
    return Chem.MolToSmiles(rdkit_mol)


def smileswriter(rdkit_mol, fname):
    
    f = open(fname, 'a')

    try:
        smiles = smilesprinter(rdkit_mol)

        line = smiles + '\n'
        f.write(line)

    except Exception as e:
        print(str(e))

        
def sdprinter(rdkit_mol, fname = ''):
    
        Chem.rdDepictor.Compute2DCoords(rdkit_mol)
        sd = Chem.MolToMolBlock(rdkit_mol)

        sd += "\n$$$$"
        return sd
 
          
def sdfilewriter(rdkit_mol, fname):
    
        f = open(fname, 'a')

        sd = sdprinter(rdkit_mol)
        f.write(sd + '\n')


# In[7]:


def initmolFromSmiles(smiles):
    """
    initialize the molecule with a structure via a smiles string.
    This allows creating all possible decorations of the scaffold structure
    provided as an argument
    """
    rdkit_mol = Chem.MolFromSmiles(smiles)
    sd = Chem.MolToMolBlock(rdkit_mol)

    sdlines = sd.split('\n')
    natoms = int(sdlines[3][0:4])
    nbonds = int(sdlines[3][4:8])
    atoms = []
    bonds = []
    
    for i in range(4,natoms+4):
        atoms.append(sdlines[i][30:34].strip())

    for i in range(4+natoms, 4+natoms+nbonds):
        line = sdlines[i]
        frm =  int(line[0:3]) - 1 # make 0 indexed
        to   = int(line[3:6]) - 1 
        order= int(line[6:9])
        chiral = int(line[9:12])
        bond = fixbond((frm, to, order, chiral))
        bonds = addbond(bonds, bond)
        
    bonds.sort()
    return (atoms, bonds)     


# In[1]:


def initmol(smiles = ''):
    
    if smiles != '':
        return initmolFromSmiles(smiles)
    else:
        # initialize with an atom
        bonds = []
        atoms = ['C']
        return (atoms, bonds)
    

def makescaffold(mol, natoms):
    """ 
    create a random scaffold of natoms carbon atoms, without rings
    """
    if natoms < 1:
        return mol
    
    atoms, bonds = mol
    num_atoms = len(atoms)
    
    # build atom-by-atom
    for bondfrom in range(num_atoms, num_atoms + natoms): # this will be the new atom number
  
        while True: # try until a bond is found somewhere
            # randomly select an atom to connect to
            bondto = random.randrange(0,bondfrom)
            if atoms[bondto] not in ringatoms:  # limit
                continue
                
            bond = fixbond((bondfrom, bondto, 1, 0))
            
            if testbond(bonds, bond)                    and valence(bonds, bondto)   < 4:
                atoms.append('C')
                bonds = addbond(bonds, bond)
                break  # found an atom to bond to
    
    return(atoms, bonds)



def gaussrandom(ave, sd):
    """
    gaussian distribution of positive integers.
    this is clunky by produces the correct distribution
    greater than 0
    """
    while True:
        result = random.gauss(ave, sd)
        if result > 0:
            return int(result)
      

def makexample(natoms,smiles=''):
       
    mol = initmol(smiles)
    natoms -= 1

    # add aromatic ring as "root" structure
    if random.random() < PROB_AROMATIC:
        mol = makearomatic(mol)
        natoms -= 6
    elif random.random() < PROB_NAPHTHYL:
        mol = makenaphthyl(mol)
        natoms -= 10

    # add some non-ring atoms
    nf = math.floor(natoms/2)
    mol = makescaffold(mol, nf)
    natoms -= nf
    
    # sometimes add another aromatic ring
    if random.random() < PROB_SECOND_AROMATIC:
        mol = makearomatic(mol)
        natoms -= 6

    # after possibly adding another aromatic ring,
    # add the rest of the atoms.
    mol = makescaffold(mol, natoms)

    # add some random ring bonds with Gaussian probability
    mol = ringbonds(mol, gaussrandom(0,1))
    atoms, bonds = mol

    #upgrade some random bonds to triple
    ntriplebonds = gaussrandom(0, 0.25)
    mol = multiplebonds(mol, ntriplebonds, 3)

    # upgrade some random bonds to double
    ndoublebonds = gaussrandom(0, len(bonds)/2)
    mol = multiplebonds(mol, ndoublebonds, 2)
    
    # randomly add some carbonyl groups for aldehydes
    # and ketones.
    if add_hetero:
        mol = makecarbonyl(mol)

        if random.random() < PROB_CARBOXYLIC:
            mol = makecarboxylic(mol)

        atoms, bonds = mol

        mol = heteroatoms(mol)

    return mol

    
    
def makemol(nmols, writer, fname, smiles='', low=12, high=25):
    seed = 198733
    random.seed(seed) # try to make it repeatable
    minimum = 1e9
    maximum = 0
    
    hash = {}

    # delete output file if exists 
    try:
        if os.path.exists(fname):
            os.remove(fname)
    except:
        pass
            
    start = time.time()
    reportinterval = 1000000
    
    for i in range(nmols):
        
        # select size of molecules
        natoms = random.randrange(low, high)
        mol = makexample(natoms, smiles)
        sd = sdfile(mol)
        rdkit_mol  = Chem.MolFromMolBlock(sd)
         
        try:
            inchi = Chem.inchi.MolToInchi(rdkit_mol, logLevel= None , treatWarningAsError=False)
        except:
            # ignore errors
            continue
        
        # keep track of molecules already generated
        if inchi in hash:
            hash[inchi] += 1
        else:
            hash[inchi] = 1
            writer(rdkit_mol, fname)
        
        # periodic statistics
        if i % reportinterval == 0 and i > 0:
            molspersecond = reportinterval/(time.time() - start)
            start = time.time()
            print('count %8d  mol/sec %4.2f'% (i, molspersecond) )

        
        if i % 10000 == 0:
            
            minimum = 100000
            maximum = -1
            
            for i in hash:
                item = hash.get(i)
                if item < minimum:
                    minimum = item
                if item > maximum:
                    maximum = item
             
            # in theory if all molecules are found N many times, the chance that
            # there are others not yet found is 1 in 2^N
            # this can be used for limited size molecules to generate complete sets
            if minimum > 6:
                break
                
            #print('hash size', len(hash), minimum, maximum)
            
    print('hash size', len(hash), 'min', minimum, 'max', maximum)        
   


def findmol(drug):
    """
    find a specific drug by limiting the generation to molecules with the
    correct number of atoms
    """
    seed = 1131317
    random.seed(seed) # try to make it repeatable
    count = 0
    
    rdkit_drug = Chem.MolFromSmiles(drug)
    drug = Chem.MolToSmiles(rdkit_drug, canonical = True)
    nheavy = rdkit_drug.GetNumAtoms()
    print(drug, 'atoms', nheavy, 'seed', seed)
    start = time.time()
    reportinterval = 1000000
    
    while True:
        count += 1
        if count % reportinterval == 0:
            molspersecond = reportinterval/(time.time() - start)
            start = time.time()
            print('count %8d  mol/sec %4.2f'% (count, molspersecond) )
            
        mol = makexample(nheavy)
        sd = sdfile(mol)
        
        rdkit_mol = Chem.MolFromMolBlock(sd)
        smiles = Chem.MolToSmiles(rdkit_mol, canonical = True)

        if smiles == drug:
            print("found", drug, ' attempt', count)
            break


# In[ ]:




