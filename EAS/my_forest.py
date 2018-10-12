 # -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 13:03:34 2018


Author: Anna Tomberg

This code is used to train and test a random forest with chosen parameters.
Example of usage:
    
    python my_forest.py train my_data.csv
    
"""



# LOAD LIBRARIES

import numpy
import pickle
import MAP_paths



def run_my_forest(X_data):  
     
    
    # LOAD DATA
    x = numpy.delete(X_data, len(MAP_paths.ALL_PROPERTIES), axis=1)
    
    print('Running RF model')
    print('*-----------------------------*')
        
    my_forest = pickle.load(open(MAP_paths.RANDOM_FOREST_MODEL, 'rb'))
    print('Loaded Random Forest Classifier from:', MAP_paths.RANDOM_FOREST_MODEL )
        
                
    m = 'Summary of Results. \n'         
    m = m + '*-----------------------------*\n'
        
    predicted_y= (my_forest.predict(x))
    class_probability = (my_forest.predict_proba(x))

    #print(str(X_data[0]))
    m = m + 'C_IDX,  PRED_ACTIVITY,   PROBAB_of_INACT, PROBAB_of_ACT \n'
    for i in range(len(predicted_y)):
        m_i = ', '.join([str(int(X_data[i][-1])),  str(predicted_y[i]), "{0:.2f}".format(class_probability[i,0]), "{0:.2f}".format(class_probability[i,1])])
        m = m + m_i + '\n'
    
    #print(m)
    
    fo = open('my_results.txt', 'w')
    fo.write(m)
    fo.close()
    print('>> Results written to my_restults.txt.')
    return predicted_y, class_probability
        


if __name__ == "__main__":
    # execute only if run as a script
    import sys
    

    filename = sys.argv[1]
    path_to_model = sys.argv[2]
    
    run_my_forest(filename, path_to_model)

    print('> DONE.')

