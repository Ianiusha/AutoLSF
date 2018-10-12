# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 10:29:53 2018

@author: RegioSQM 
Alternations by A. Tomberg
"""

from MAP_gen_g09 import MAP_generate_g09
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
import re
from openbabel import OBConversion, OBMol
import MAP_paths



def get_num_hydrogens(atom):
    return len([True for a in atom.GetNeighbors() if a.GetSymbol() == 'H'])


def check_if_aromaticCH(my_mol):
    found_aromatic_CH = False
    m = Chem.AddHs(my_mol)
    arom_carbons =  m.GetAromaticAtoms()
    for atom in arom_carbons:
        if atom.GetSymbol() == 'C':
            
            if get_num_hydrogens(atom) == 1:
                found_aromatic_CH = True
	
    return found_aromatic_CH


def convert_g09_to_SDF(filename):
    
    
    obConversion = OBConversion()
    obConversion.SetInAndOutFormats("g09", "sdf")
    #obConversion.AddOption("-h")

    mol = OBMol()
    obConversion.ReadFile(mol, filename)  

    new_file = 'original.sdf'
    obConversion.WriteFile(mol, new_file)
    
    return new_file




def generate_conformations(m, n):
    
    mol = Chem.AddHs(m)
    ids=ids = AllChem.EmbedMultipleConfs(mol,numConfs=n,useExpTorsionAnglePrefs=True,useBasicKnowledge=True)
    #ids=ids = AllChem.EmbedMultipleConfs(mol,numConfs=n)
    results ={}
    
    for i in ids:   
        
        try:
            
            ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=i)
            ff.Initialize()
            ff.CalcEnergy()
            if MAP_paths.MINI_Iterations > 0:
                AllChem.MMFFOptimizeMolecule(mol, confId=i)
                results[i] = ff.CalcEnergy()
                
        except:
            
            ff = AllChem.UFFGetMoleculeForceField(mol, AllChem.UFFGetMoleculeProperties(mol), confId=i)
            ff.Initialize()
            ff.CalcEnergy()
            if MAP_paths.MINI_Iterations > 0:
                AllChem.UFFOptimizeMolecule(mol, confId=i)
                results[i] = ff.CalcEnergy()

    return mol, results




def generate_lowestE_conformer(my_smiles):
     
    m = Chem.AddHs(Chem.MolFromSmiles(my_smiles))

    # DECIDE how many conformers are to be computed
    rot_bond = rdMolDescriptors.CalcNumRotatableBonds(m)
    n_confs = min(1 + 3*rot_bond,MAP_paths.MAX_CONF)
    
    # GENERATE conformers
    my_mol, confID_energies = generate_conformations(m, n_confs)
    
    # find lowest energy ID
    min_conformer_id = min(confID_energies, key=confID_energies.get)
        
    return my_mol, min_conformer_id


    

def get_energy(my_name):
    
    H_energy = 0
    out_block = ''
    solvationE = 0
    output_name = my_name + '_PM3.log'
    print(output_name)
    total_E = 'NA'

    fo = open(output_name, 'r')
    for line in fo.readlines():
        if re.search('\\\\', line):
            out_block = out_block + (line.rstrip()).replace(" ", "")
        # look for solvation energy:
        elif re.search('DeltaG \(solv\)', line):
            solvationE = float(line.split()[-1])
    fo.close()

    if not (out_block == ''):
        t = (out_block.split('\\HF='))
        if len(t) == 2:
            a = (t[1].split('\\'))[0]
            # look for a second energy value
        elif len(t) == 3:
            a = (t[2].split('\\'))[0]

        H_energy = float(a)

        if (H_energy != 0) and (solvationE != 0):
            total_E = (H_energy * MAP_paths.CONVERT_AU2KCAL_MOL + solvationE)
            # print  (gname, H_energy , solvationE, total_E)

    return total_E


def generate_charged_smiles(my_molecule):
    """ 
    Here we will...
        - generate a list of protonates states for each aromatic carbon with a H
        - return a dictionary of compounds, with corresponding protonated states
    """    

    name_list = []
    smiles_list = []
    atom_list = []

   
    arom_carbons_idxs = [a.GetIdx() for a in my_molecule.GetAromaticAtoms()]

    Chem.Kekulize(my_molecule,clearAromaticFlags=True)

    counter = 0
    for idx in arom_carbons_idxs:
        
        atom = my_molecule.GetAtomWithIdx(idx)
        
        if atom.GetSymbol() == 'C':

            
            if  get_num_hydrogens(atom) == 1:
                
                
                counter = counter + 1
            
                mw_i = Chem.RWMol(my_molecule)
                carbon = mw_i.GetAtomWithIdx(idx)
                                
                idx_H = mw_i.AddAtom(Chem.Atom(1))
                mw_i.AddBond(idx_H, idx, Chem.BondType.SINGLE)
                
                new_hydrogen = mw_i.GetAtomWithIdx(idx_H)
                new_hydrogen.UpdatePropertyCache()
                

                for bi in carbon.GetBonds():
                    if bi.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        
                        bi.SetBondType(Chem.rdchem.BondType.SINGLE)
                        bonded_to = bi.GetOtherAtom(carbon)
                        bonded_to.SetFormalCharge(1)
                        bonded_to.UpdatePropertyCache()
                
                carbon.UpdatePropertyCache()        
 
                #img = Draw.MolsToGridImage([mw_i],legends=[ "protonated on C"+str(idx)], kekulize=False)
                #img.save(str(idx)+'.png')
                
                Chem.SanitizeMol(mw_i)
                smiles = Chem.MolToSmiles(mw_i)
                   
        
                name_list.append("m_" + str(counter))
                smiles_list.append(smiles)
                atom_list.append(atom.GetIdx())     
                

    return name_list, smiles_list, atom_list



def write_conformer_to_sdf(mol,conformer_id, filename):
    
    tm = Chem.Mol(mol, False, conformer_id)
    w = Chem.SDWriter(filename)
    w.write(tm)
    w.flush()
    w.close()
    
    return filename



def generate_conformations_files(smiles, name, charge, max_conf, option='PM3'):
    
    #csmile, cname, ccharge, 20

    "Use only lowest energy conformer per protonated state. "
    

    my_mol, conformer_ID = generate_lowestE_conformer(smiles)
    fname = write_conformer_to_sdf(my_mol, conformer_ID, name + '.sdf')

    if option == 'PM3':
        # Create a gaussian input for PM3 calculation
        g09_input = MAP_generate_g09(fname, charge, 'PM3')

        " RUN THE JOBS CREATED:"  
        #submit_g09_job(g09_input)
        
    elif option == 'OPT':

        # Create a gaussian input for PM3 calculation
        g09_input = MAP_generate_g09(fname, charge, 'OPT')

        " RUN THE JOBS CREATED:"  
        #submit_g09_job(g09_input)
    return 
