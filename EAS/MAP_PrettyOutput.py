# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 15:00:52 2018

@author: kntb373
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 10:29:53 2018

@author: RegioSQM 
Alternations by A. Tomberg
"""

import sys
import rdkit.Chem
import rdkit.Chem.AllChem
from rdkit.Chem import Draw




def generate_color_dict(sure, unsure, very_unsure):
    colors_dict = {}
    
    for a in very_unsure:
        colors_dict[a] = (1,0,0)    
    for a in unsure:
        colors_dict[a] = (1,1,0)
    for a in sure:
        colors_dict[a] = (0,1,0)
        
    return colors_dict

def prettyOutput(parent_sdf, reactivity_dictionary):


    suppl = rdkit.Chem.SDMolSupplier(parent_sdf)
    mol = suppl[0]
    mol = rdkit.Chem.RemoveHs(mol)
    rdkit.Chem.Kekulize(mol)
    tmp = rdkit.Chem.AllChem.Compute2DCoords(mol)  
    
    sure, unsure, very_unsure = prep_atom_list(reactivity_dictionary)
    
    colors = generate_color_dict(sure, unsure, very_unsure)
    merged_atoms = sure+unsure+very_unsure

    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(500,500)
    op = drawer.drawOptions()

    for highlighted_atom in reactivity_dictionary.keys():
       op.atomLabels[highlighted_atom]= (mol.GetAtomWithIdx(highlighted_atom)).GetSymbol() + str((highlighted_atom+1))
    
    drawer.DrawMolecule(mol,highlightAtoms= merged_atoms, highlightAtomColors=colors)

       
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg2 = svg.replace('svg:', '')


    return svg2



def prep_atom_list(reactivity_dictionary):
 
    max_probab = max(list(reactivity_dictionary.values()))
    #print(max_probab)
    very_uncertain = []
    uncertain = []
    sure = []
    
    for idx, probab in reactivity_dictionary.items():
        
        if probab >= 0.5 and probab < 0.70:
                uncertain.append(idx)
        elif probab >= 0.70:
                sure.append(idx)

    
    if max_probab < 0.50:
        # all carbons were found to be inactive, 
        # make the highest probab carbon active
        idx = (list(reactivity_dictionary.keys())[list(reactivity_dictionary.values()).index(max_probab)])
        very_uncertain.append(idx)        
                
    #predicted = (sure, uncertain)            
    #print(predicted)            
    return sure, uncertain, very_uncertain


def svg_with_atom_numbers(my_sdf, keep_hydrogens=False, special_atoms = None):

    
    if keep_hydrogens: 
        suppl = rdkit.Chem.SDMolSupplier(my_sdf, removeHs=False)
        mol = suppl[0]
    else:
        suppl = rdkit.Chem.SDMolSupplier(my_sdf)
        mol = suppl[0]
        mol = rdkit.Chem.RemoveHs(mol)
        
    
    rdkit.Chem.Kekulize(mol)
    tmp = rdkit.Chem.AllChem.Compute2DCoords(mol)    

    dr = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DSVG(500,500)
    op = dr.drawOptions()
    
    if special_atoms == None:
        special_atoms = mol.GetAtoms()
        
    for special_atom in special_atoms:
       op.atomLabels[special_atom]= (mol.GetAtomWithIdx(special_atom)).GetSymbol() + str((special_atom+1))


    dr.DrawMolecule(mol)
    dr.FinishDrawing()
    svg = dr.GetDrawingText()
    svg2 = svg.replace('svg:', '')
    
    # format for XHTML:
    my_svg = svg2.split('\n')[7:]    
    my_svg = ''.join(my_svg)  
    #print(my_svg)
 
    return my_svg

def basic_svg_from_smi(my_smiles):

    mol = rdkit.Chem.MolFromSmiles(my_smiles)
    
    rdkit.Chem.Kekulize(mol)
    rdkit.Chem.AllChem.Compute2DCoords(mol)
    

    dr = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DSVG(500,500)
    dr.DrawMolecule(mol)
    dr.FinishDrawing()
    svg = dr.GetDrawingText()
    svg2 = svg.replace('svg:', '')
    
    # format for XHTML:
    my_svg = svg2.split('\n')[7:]    
    my_svg = ''.join(my_svg)  
    #print(my_svg)
 
    return my_svg


if __name__ == "__main__":
    # execute only if run as a script

    option = sys.argv[1]
    filename = sys.argv[2]

    if option == 'enumerated':
        svg_with_atom_numbers(filename)



    print('> DONE.')
