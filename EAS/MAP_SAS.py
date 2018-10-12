# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 10:29:53 2018

@author: kntb373
"""



import MAP_gen_g09
import  MAP_paths 
import subprocess



def reformat_gaussian_xyz(xyz_block):
    # the xyz block that comes out from the gaussian log parser is 
    # in the wrong format for GEPOL; reformat and keep as list to add the 
    # radii later
  
    new_xyz_list = []
    #vdw_radii = get_radii_dictionary()
    vdw_radii = MAP_paths.VDW_DICTIONARY
    
    mol_list = xyz_block.split('\n')
    for m in mol_list[:-1]:
        
        #print (m)
        m = m.split()
        a = (m[0])
        x = float(m[1])
        y = float(m[2])
        z = float(m[3])
        
        if a in vdw_radii:
            r = vdw_radii[a]
        else:
            print('>> WARNING: the atom ' + a \
                  +' is not in VDW_DICTIONARY. It got assigned a default value of '\
                  + str(MAP_paths.DEFAULT_VDW_RADIUS))
            
            r = float(MAP_paths.DEFAULT_VDW_RADIUS)
        
        
        #'{: f}; {: f}'.format(3.14, -3.14)
        m = ' {: .5f}  {: .5f}  {: .5f}  {: .5f}     {}'.format(x,y,z,r,a)
        new_xyz_list.append(m)
    #print (new_xyz_list)
    return new_xyz_list


"""
def get_radii_dictionary():
    # we need to form a dictionary with atoms and the corresponding 
    # van der Waals radii. Taking from/gepol93/VdWradii.lib
    
    vdw_dictionary = {}
    
    try:
            fo = open(MAP_paths.VdW_LIBRARY , 'r')
            data = fo.readlines()
            fo.close()
            
            for l in data:
                #print(l)
                if not(l.startswith('*') or l in ['\n', '\r\n']):

                    atom = l.split()[0]
                    radius = l.split()[-1]
                    vdw_dictionary[atom] = float(radius)
                    #print(vdw_dictionary)

    except:
            print (">> ERROR: cannot get radii from library: " +MAP_paths.VdW_LIBRARY)
            return 
    

    #print(vdw_dictionary)
    return vdw_dictionary
"""

def make_GEPOL_inputs(base_name, xyz_block, solvent_radius = None):
    # make 2 files required for GEPOL
    # (1) control file : keywords for SAS computation (*.INP)
    # (2) coordinates file : xyz and radii list (*.XYZR)
    
    if solvent_radius == None:
        solvent_radius = '1.4'
    
    control_file = base_name + '.INP'
    coordinates_file = base_name + '.XYZR'

    
    # (1) control file : keywords for SAS computation (*.INP)   
    my_input =  MAP_paths.SAS_INP.substitute(NAME=base_name, SOLVENT_RAD = solvent_radius)
    fo = open(control_file, 'w')
    fo.write(my_input)
    fo.close()
    
    print('>> ' + control_file)

    # (2) coordinates file : xyz and radii list (*.XYZR)
    my_input = MAP_paths.SAS_XYZR.substitute(NAME=base_name, ATOM_NUMBER=len(xyz_block))
    for i in xyz_block:
        my_input = my_input + i + '\n'
        
    fo = open(coordinates_file, 'w')
    fo.write(my_input)
    fo.close()
    
    print('>> ' + coordinates_file)
    
    return control_file, coordinates_file


def get_SAS(filename, molecule_container):
    """
    Will set the following atomic descriptors for all non-Hydrogens:
        SAS = 
    """ 
    total_atoms = molecule_container.GetNumAtoms()


    fo = open( 'mySAS.SPH', 'r')
    data = fo.readlines()
    fo.close()

    try:
        for i in range (total_atoms):
            l = data[i+3]
            if not( l in ['\n', '\r\n']):
                temp = (l).split()
                #print(temp)
               
                current_atom = molecule_container.GetAtomWithIdx(i)
                #print(current_atom.GetSymbol())
    
                if current_atom.GetSymbol() == temp[4]:
                    # Index corresponds to atom 
                    current_atom.SetProp('SAS', temp[-1])
    except:
        print (">> ERROR: cannot read atomic SAS.")
        
    
    return molecule_container




def MAP_computeSAS(filename, option=1.4):
    
    """
        Here we will ...
            - check what type of input file is given
            - will use only the opted geometry for SAS
    """    
    solvent_radius = option
    
    #1) figure out how many molecules are initially in workspace
    #       ... 1 if gaussian log or sdf
    #       ... multiple if smiles file
       
        
    print(">> Running SAS from GEPOL with gaussian input.")
        
    """
        - For only the last geometry in the log file...
        ... prepare GEPOL input files
        ... add their names to the "to_run" list
        ... run GEPOL
        ... collect results and parse into the feature file
        
        * WORKING WIHT ONLY 1 CONFORMER *
    """
        
        #1) extract xyz coords and atom list from gaussian log
    charge,multiplicity,coordinates = MAP_gen_g09.parse_gaussianLOG(filename)
    new_coordinates = reformat_gaussian_xyz(coordinates)
        
        #2) create GEPOL input files
    INP, XYZR = make_GEPOL_inputs('mySAS' , new_coordinates, solvent_radius)
        
        #3) run GEPOL93 SAS job
    subprocess.call('sas.out < '+INP + ' > ' + filename.replace('.log', '.gepol93'), shell=True)        
                
                
        
    return
    
    





 
