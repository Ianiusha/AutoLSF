# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 10:29:53 2018

@author: kntb373

Routine:
    1.	Convert original.sdf to molecule object
    2.	Extract values of different descriptors.
    3.	Compile data into a summary file
    4.	Prepare input for random forest
    5.	Run random forest and output prediction


"""

import numpy as np
import re
from rdkit import Chem

from rRSQM_support import get_energy, get_num_hydrogens
import MAP_SAS 
import MAP_paths
from my_forest import run_my_forest

from MAP_PrettyOutput import prettyOutput


def get_DDEC_Qs(molecule_container):
    """
    Will set the following atomic descriptors for all non-Hydrogens:
        Q = 
        Bond_orders_sum = 
    """
    total_atoms = molecule_container.GetNumAtoms()
    # print(total_atoms)

    fo = open('DDEC6_even_tempered_net_atomic_charges.xyz', 'r')
    data = fo.readlines()
    fo.close()

    try:
        for i in range(total_atoms):

            temp = (data[i + 2]).split()
            # print(temp)

            current_atom = molecule_container.GetAtomWithIdx(i)
            # print(current_atom.GetSymbol())

            if current_atom.GetSymbol() == temp[0]:
                # Index corresponds to atom 
                current_atom.SetProp('Q', temp[-1])
    except:
        print(">> ERROR: cannot read atomic charges.")

    print("Extracted DDEC partial charges.")
    return molecule_container


def get_DDEC_BOs(molecule_container):
    """
        Will set the following atomic descriptors for all aromatic carbons with Hs:
        CH_Bond_order = 
    """
    total_atoms = molecule_container.GetNumAtoms()
    fo = open('DDEC6_even_tempered_bond_orders.xyz', 'r')
    data = fo.readlines()
    fo.close()

    try:
        
        all_aromatic_carbons = [x.GetIdx() for x in molecule_container.GetAromaticAtoms()]
        
        for i in range(total_atoms):

            temp = (data[i + 2]).split()
            # print(i, temp)

            current_atom = molecule_container.GetAtomWithIdx(i)
            # print(current_atom.GetSymbol())

            if current_atom.GetSymbol() == temp[0]:
                # Index corresponds to atom 
                current_atom.SetProp('Bond_orders_sum', temp[-1])

                if current_atom.GetIdx() in all_aromatic_carbons:

                    # get hydrogen that is bonded to this aromatic carbon
                    my_hydrogens = []
                    for n in current_atom.GetNeighbors():
                        # print(n.GetSymbol())
                        if n.GetSymbol() == 'H':
                            my_hydrogens.append(str(n.GetIdx() + 1))
                    # print(my_hydrogens)

                    if len(my_hydrogens) > 0:
                        # this aromatic carbon has a hydrogen 

                        counter = 0
                        for line in data:
                            if re.search(' Printing BOs for ATOM # (\s)* ' + str(i + 1), line):
                                # print(line)
                                break
                            counter = counter + 1

                        # print( current_atom.GetSymbol(), str(i+1), counter, my_hydrogens)
                        line = data[counter]
                        c = 0
                        while not (line.startswith(' ======================')):
                            # print(line)
                            if line.startswith(' Bonded to the (  0,   0,   0) translated image of atom number'):
                                temp = line.split()
                                if temp[12] in my_hydrogens:
                                    # print(temp[20])
                                    current_atom.SetProp('CH_Bond_order', temp[20])
                            c = c + 1
                            line = data[counter + c]

    except:
        print(">> ERROR: cannot read ond orders.")

    print("Extracted DDEC bond orders.")
    return molecule_container


def get_FMO_fukuis(molecule_container):
    """
    Will set the following atomic descriptors for all non-Hydrogens:
        FUKUI = 
    """
    total_atoms = molecule_container.GetNumAtoms()

    fo = open('overlap_populations.xyz', 'r')
    data = fo.readlines()
    fo.close()

    # create a 2D matrix for the overlap between atoms
    overlap_matrix = np.zeros((total_atoms, total_atoms), dtype=float)
    # print('OVERLAP')
    total_integrals = int(data[0])
    for line in data[3:(total_integrals + 3)]:
        # print('>', line)
        t = line.split()
        atom1 = int(t[0]) - 1
        atom2 = int(t[1]) - 1
        overlap = float(t[-1])
        overlap_matrix[atom1, atom2] = overlap
        # print(atom1, atom2, overlap)

    # Get the HOMO densities from corresponding Gaussian log file
    fo = open('parent_OPT.log', 'r')
    data = fo.readlines()
    fo.close()

    # Determine the # of the HOMO orbital
    for i in range(len(data)):
        line = data[i]
        m = re.search('(\d+) alpha electrons(\s)+(\d+) beta electrons', line)
        if m:
            # print('>', line, m.groups())
            if m.groups()[0] == m.groups()[2]:
                HOMO_nbr = int(m.groups()[0])
                # print('HOMO #:', HOMO_nbr)

            else:
                print('>> WARNING: found a different # of alpha and beta electrons!')

        elif line.startswith(' Alpha  occ. eigenvalues --'):
            HOMO_energy = line.split()[-1]
            # print('HOMO E:', HOMO_energy)

        elif re.search('Molecular Orbital Coefficients: ', line):
            break

    # initiate HOMO_density array
    HOMO_density = np.zeros((total_atoms))
    # homo_row = ceil(HOMO_nbr/5) # This give section where the homo coeeficients start

    p = re.compile('(Eigenvalues --)(.)+ ' + HOMO_energy)
    for i in range(len(data)):
        m = p.search(data[i])
        if m:
            l = (data[i]).split()
            index = 7 - l.index(HOMO_energy)
            line_number = i
            # print('>',i, l) # i = line where the HOMO starts
            # print('>',index, data[i-2])
            # sanity check:
            if not ((data[i - 2]).split()[-index] == str(HOMO_nbr)):
                print('>> WARNING: HOMO density confusion!')
                return -1

            break

    # print(line_number)
    counter = line_number
    for line in data[(counter + 1):]:
        counter = counter + 1
        l = line.split()
        # print(l)
        # print(data[counter+1])
        # first line for an atom
        if len(l) == 9:
            current_atom = int(l[1]) - 1
            # print(current_atom)
            HOMO_density[current_atom] = 0

            # print(current_atom, l[0], float(l[-index]))
        HOMO_density[current_atom] = HOMO_density[current_atom] + float(l[-index]) * float(l[-index])

        if (data[counter + 3]).startswith('     Eigenvalues --'):
            break

    # print('FUKUI')
    for i in range(total_atoms):
        current_atom = molecule_container.GetAtomWithIdx(i)

        fukui = HOMO_density[i]
        for j in range(total_atoms):
            fukui = fukui + HOMO_density[i] * HOMO_density[j] * overlap_matrix[i, j]
            # print(HOMO_density[i], HOMO_density[j], overlap_matrix[i,j])

        fukui = "{:6.6f}".format(fukui)

        current_atom.SetProp('Fukui', fukui)
        # print(i+1, current_atom.GetSymbol(),fukui)

    print("Extracted Condensed fukui coeffs.")
    return molecule_container


def get_PM3_energies(my_molecule):
 
    # Read protonates conformations list 
    f = open('m.files', 'r')

    # skip header
    f.readline()
    p = f.readline()
    my_compound_dictionary = {}
    problematic_centers = []
    # read m.files
    for line in f:

        line = line.split(", ")
        name = line[0]
        reaction_center = int(line[2])

        heat = get_energy(name)

        if heat == 'NA':
            # something went wrong with the PM3 calculation, cannot include 
            # this atom in the dictionary
            print('>> WARNING: cannot get energy from', name, '!')
            problematic_centers.append(reaction_center)
            

        else:
            my_compound_dictionary[heat] = [name, reaction_center]

    f.close()

    all_heats = np.array(list(my_compound_dictionary.keys()))
    minimum_heat = np.min(all_heats)
    
    # Assign energy to each C-H pair in my_molecule
    # for each reaction center, find all heats, and take the minimum one, assign to that center

    for k,v in my_compound_dictionary.items():
        
        #print (k, v)
        current_heat = k - minimum_heat
        #n = v[0]
        rxn_center = v[1]

        current_atom = my_molecule.GetAtomWithIdx(int(rxn_center))
        current_atom.SetProp('CH_rSQM_E', str(("%.4f" % current_heat)))
        # print (a, current_heats)
    
    for c in problematic_centers:
        current_atom = my_molecule.GetAtomWithIdx(int(c))
        current_atom.SetProp('CH_rSQM_E', str(("%.4f" % 30)))
        

    return my_molecule


def make_feature_file(molecule_container):
    
    # Makes a csv file for this molecule where each atom is an entry
    # if file already exists, adds to what's there
    feature_file = 'my_features.csv'

    my_features = 'Atom#, Atom, ' +  (', '.join(MAP_paths.ALL_PROPERTIES)) + '\n'
    
    molecule_container = Chem.RemoveHs(molecule_container)
    total_atoms = molecule_container.GetNumAtoms()

    # PREPARE DESCIPTORS CONTAINER:
    my_descriptors = []
    for atom in molecule_container.GetAtoms():
        d = {'IDX': atom.GetIdx(), 'Symbol': atom.GetSymbol()}

        for prop in MAP_paths.ALL_PROPERTIES:
            d[prop] = 'NA'

        my_descriptors.append(d)

    print("Looking for equivalent atoms...")
    Chem.AssignStereochemistry(molecule_container, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    #equivalencies = [atom.GetProp('_CIPRank') for atom in m.GetAtoms()]
    #print (equivalencies)

    found_equivalent = []

    for atom in molecule_container.GetAtoms():
        #print (atom.GetIdx(), atom.GetProp('_CIPRank'))
        if atom.GetSymbol() == 'C':
            for a in molecule_container.GetAtoms():
                if atom.GetIdx() != a.GetIdx():
                    if atom.GetProp('_CIPRank') == a.GetProp('_CIPRank'):
                        # found equivalent atoms, save index
                        found_equivalent.append([atom.GetIdx(), a.GetIdx()])

    for a in found_equivalent:
        a.sort()

    found_equivalent = set(map(tuple, found_equivalent))

    # print("EQ:", found_equivalent )
    # If CIPRank is same - atoms are equivalent! Need to average their descriptors
    # """
    for i, j in found_equivalent:

        atom1 = molecule_container.GetAtomWithIdx(i)
        atom2 = molecule_container.GetAtomWithIdx(j)

        for prop in MAP_paths.ALL_PROPERTIES:
            if atom1.HasProp(prop) and atom2.HasProp(prop):
                temp1 = float(atom1.GetProp(prop))
                temp2 = float(atom2.GetProp(prop))
                temp = (temp1 + temp2) / 2.0
                temp = "{:6.6f}".format(temp)

                atom1.SetProp(prop, temp)
                atom2.SetProp(prop, temp)
                print("AVERAGED", prop, "for", i, j)
    # """

    for i in range(total_atoms):
        current_atom = molecule_container.GetAtomWithIdx(i)
        my_features = my_features + str(current_atom.GetIdx()) + ', ' + current_atom.GetSymbol()

        for prop in MAP_paths.ALL_PROPERTIES:
            if current_atom.HasProp(prop):
                temp = current_atom.GetProp(prop)
            else:
                temp = 'NA'
            my_features = my_features + ', ' + temp

        my_features = my_features + ' \n'

    fw = open(feature_file, 'w')
    fw.write(my_features)
    fw.close()
    # print(my_features)

    print("Wrote into a features file.")
    return feature_file


def extractDescriptors(my_sdf):
    """
    -extr = extract_descriptors
    
    Here we will ...
        - find the gaussian output or other output file corresponding to
          the option selected
        - parse the file and extract the information required
        - add the information to a feature file for the molecule
    """

    # 1. convert sdf file into molecule object

    my_molecule = (Chem.SDMolSupplier(my_sdf))[0]
    my_molecule = Chem.AddHs(my_molecule)
    # 2. extract necessary data and add properties to molecule object

    my_molecule = get_DDEC_Qs(my_molecule)
    my_molecule = get_DDEC_BOs(my_molecule)
    my_molecule = get_PM3_energies(my_molecule)
    my_molecule = MAP_SAS.get_SAS(my_sdf, my_molecule)
    my_molecule = get_FMO_fukuis(my_molecule)


    # 3. compile data together
    make_feature_file(my_molecule)
    
    my_data = compile_data_dictionary(my_molecule)
    
    return my_data


def compile_data_dictionary(my_molecule):
    

    #my_molecule = Chem.RemoveHs(my_molecule)
    descriptors_data = { 'idx':[]}
    
    for d in MAP_paths.ALL_PROPERTIES:
        descriptors_data[d]=[]
    
    arom_carbons =  my_molecule.GetAromaticAtoms()
    
    #print(arom_carbons)
    l = 0
    for atom in arom_carbons:
        if atom.GetSymbol() == 'C':
            
            if  get_num_hydrogens(atom) == 1:
                #print(atom.GetIdx())
                descriptors_data['idx'].append(atom.GetIdx())
                descriptors_data['Q'].append(atom.GetProp('Q'))
                descriptors_data['CH_Bond_order'].append(atom.GetProp('CH_Bond_order'))
                descriptors_data['Bond_orders_sum'].append(atom.GetProp('Bond_orders_sum'))
                descriptors_data['SAS'].append(atom.GetProp('SAS'))
                descriptors_data['Fukui'].append(atom.GetProp('Fukui'))
                descriptors_data['CH_rSQM_E'].append(atom.GetProp('CH_rSQM_E'))
                l= l+1
    # add descriptors columns to X in same order as in all_descriptors
    #labels = numpy.array(descriptors_data['labels'])
    
    X = np.zeros((len(MAP_paths.ALL_PROPERTIES)+1 , l))
    
    
    for prop in MAP_paths.ALL_PROPERTIES:
        #print (prop,  MAP_paths.ALL_PROPERTIES.index(prop))
        X[MAP_paths.ALL_PROPERTIES.index(prop):] = descriptors_data[prop]
        
    X [(len(MAP_paths.ALL_PROPERTIES)):] = descriptors_data['idx']
    X = X.transpose()

    return X #, labels


if __name__ == "__main__":

    MAP_SAS.MAP_computeSAS('parent_OPT.log')   
    my_data = extractDescriptors('parent.sdf')
    predicted_y, class_probability = run_my_forest(my_data)
    #Draw resulting molecule
    reactivity_dictionary = {}
            
    for i in range(len(predicted_y)):
        probab_of_active = (class_probability[i][0])
        reactivity_dictionary[i] = probab_of_active
                
    # make color coded svg
    labeled_svg = prettyOutput('parent.sdf', reactivity_dictionary)
    with open('labeled_molecule.html', 'w') as f:
        f.write(labeled_svg)
    print(">> Made 2D representation of final molecule in labeled_molecule.html.")
    
    