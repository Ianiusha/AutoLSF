# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 10:11:37 2018

@author: Anna Tomberg


"""

def main_01(input_mol):
    """
    Routine:

1) convert input into rdkit molecule
2) check if arom. CH exist
3) create sdf for paretn molecule and protonated conformers
4) run gaussian jobs to compute QM descriptors

    """
    
    import subprocess
    from rdkit import Chem
    import rRSQM_support

    
    #1. figure out if smiles or sdf
    converted_to_rdkit_mol = False
    try:
        my_molecule = (Chem.SDMolSupplier(input_mol))[0]
        converted_to_rdkit_mol = True
    except:
        converted_to_rdkit_mol = False
    
    if not(converted_to_rdkit_mol):
        try:
            my_molecule = Chem.MolFromSmiles(input_mol)  

        except subprocess.CalledProcessError as exc:
            print(">> ERROR: cannot read input molecule.", exc.returncode, exc.output)
            return

    # CONVERT BACK TO SMILES AND THEN SDF TO HAVE SAME NUMBERING 
    my_smiles = Chem.MolToSmiles(my_molecule,isomericSmiles=True)
    my_molecule = Chem.AddHs(Chem.MolFromSmiles(my_smiles))

    if rRSQM_support.check_if_aromaticCH(my_molecule):
        

        # 2. run protonation protocol
        cnames, csmiles, catoms  = rRSQM_support.generate_charged_smiles(my_molecule)
        print(Chem.MolToMolBlock(Chem.AddHs(Chem.MolFromSmiles(my_smiles))))
        
        if len(csmiles) < 1:
                print("ERROR: protonation protocol went wrong. Exit.")
                return
            
        pcharge = Chem.GetFormalCharge(my_molecule)
        # 3. gaussian opt job and DDEC
        rRSQM_support.generate_conformations_files(my_smiles, 'parent', pcharge, 20, option='OPT')
        
        m = "name, SMILES, reaction_center \n"
        m = m + ( ", ".join(['m', my_smiles]) +" , charge={}".format(str(pcharge)) + '\n')



        ccharge = pcharge + 1
		
        for cname, csmile, catom in zip(cnames, csmiles, catoms):
		
                # Do conformational search on each smiles structure and save it in SDF format
                rRSQM_support.generate_conformations_files(csmile, cname, ccharge, 20)
                m = m + ( ", ".join([cname, csmile, str(catom)]) +" , charge={}".format(str(ccharge)) + '\n')
                
                fw = open('m.files', 'w')
                fw.write(m)
                fw.close()
                

    else:
        print('>> There are no aromatic C-H found in this molecule. Abort.')    
        return



def main_02():

    import MAP_SAS 
    from EAS_RF import extractDescriptors
    from my_forest import run_my_forest
    from MAP_PrettyOutput import prettyOutput

    MAP_SAS.MAP_computeSAS('parent_OPT.log')   
    my_data = extractDescriptors('parent.sdf')
    run_my_forest(my_data)

    fo = open('my_results.txt', 'r')
    txt = fo.readlines()
    fo.close()
            
    reactivity_dictionary = {}
            
    for line in txt[3:]:
           line = (line.split(','))
           idx = int(line[0])
           probab_of_active = float(line[3].strip())
                
           reactivity_dictionary[idx] = probab_of_active
    # make color coded svg
    labeled_svg = prettyOutput('parent.sdf', reactivity_dictionary)
    with open('labeled_molecule.html', 'w') as f:
        f.write(labeled_svg)
    print(">> Made 2D representation of final molecule in labeled_molecule.html.")
    
    
    return




if __name__ == "__main__":
    import sys

    if sys.argv[1] == 'compute':
        # argument can be either an sdf, or a smiles string with ONE molecule
        main_01(sys.argv[2])

    elif sys.argv[1] == 'extract':
        main_02()