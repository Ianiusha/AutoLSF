# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 12:18:21 2018

Author: Anna Tomberg

This code is used to format the input for machine learning models.
It takes in a csv file with the descriptors of carbons. Format of csv:
    comp1, 1, 0, 0.005984, 0.9121, 4.094493, 17.72268, 0.010775, 9.6843, 
    name, experimental_label, Q, CHBO, SBO, SAS, FC, rRSQM,
    
Labels are 0 = inactive in EAS ; 1 = active in EAS
The data set should have no missing values.
"""



# LOAD LIBRARIES
import csv, numpy


###############################################################################

qm_descriptors =  ['Q', 'CHBO', 'SBO', 'SAS',  'FC', 'rRSQM']

###############################################################################


def load_data (data_file, descriptors):

    my_descriptors = qm_descriptors
    print(len(my_descriptors))
    
    descriptors_data = {'labels':[] , 'idx':[]}
    
    
    for i in my_descriptors:
        descriptors_data[i]=[]
    
    with open(data_file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        
        l = 0
        
        for row in readCSV:

            descriptors_data['idx'].append(l)
            
            descriptors_data['Q'].append(float(row[2]))
            descriptors_data['CHBO'].append(float(row[3]))
            descriptors_data['SBO'].append(float(row[4]))
            descriptors_data['SAS'].append(float(row[5]))
            descriptors_data['FC'].append(float(row[6]))
            descriptors_data['rRSQM'].append(float(row[7]))
            
            descriptors_data['labels'].append(float(row[8]))
      
            l = l+1
        
    #arr9,8, -0.047535, 0.891200, 4.004582, 11.573105, 0.103566, 1.211750, 1.000000 
    # add descriptors columns to X in same order as in all_descriptors
    labels = numpy.array(descriptors_data['labels'])
    
    print(labels)
    
    X = numpy.zeros((len(my_descriptors)+1 , len(labels)))
    
    
    for i in my_descriptors:
        #print (i, my_descriptors.index(i))
        X[my_descriptors.index(i):] = descriptors_data[i]
        
    X [(len(my_descriptors)):] = descriptors_data['idx']
    X = X.transpose()
    #print(X[0])

        
    print ('active LABELS: %.2f' % (100*sum(labels)/len(labels)))
    print ('Shape of X data: ', (X.shape))
    print ('Shape of Y data: ', (labels.shape))
    
    return X, labels



def do_down_sampling(x,y):
    
    #https://chrisalbon.com/machine_learning/preprocessing_structured_data/handling_imbalanced_classes_with_downsampling/

    # Indicies of each class' observations
    negatives = numpy.where(y == 0)[0]
    positives = numpy.where(y == 1)[0]

   
    # Number of observations in each class
    n_negatives = len(negatives)
    n_positives = len(positives)
    
                    
    print('>> Found', n_positives, 'positives and ', n_negatives, 'negatives.')
    print('>> Reducing data size by random sampling ...')
    
    # For every observation of class 0, randomly sample from class 1 without replacement
    negatives_downsampled = numpy.random.choice(negatives, size=n_positives, replace=False)
    
    print(len(negatives_downsampled))
    
    # Join together class 1's target vector with the downsampled class 0's target vector
    y_resampled = numpy.concatenate((y[positives], y[negatives_downsampled]))
    x_resampled = numpy.concatenate((x[positives], x[negatives_downsampled]))
    
    print('>> New size of data: ', len(y_resampled))
    
    return x_resampled, y_resampled



def split_negatives_positives(x,label):
    
    
    negatives = numpy.where(label == 0)[0]
    positives = numpy.where(label == 1)[0]
                
    positives_x = numpy.zeros(len(positives)) 
    j = 0
       
    for index in positives:
            positives_x[j] = x[index]
            j = j + 1
            
    negatives_x = numpy.zeros(len(negatives)) 
    j = 0
       
    for index in negatives:
            negatives_x[j] = x[index]
            j = j + 1
    
        
    return positives_x, negatives_x



def split_data (X, Y, partition):
    
    from sklearn.model_selection import train_test_split
    
    x1, x2, y1, y2 = train_test_split(X, Y, test_size = partition)
    
    print ('SET 1: \n\t Num.of data points: '+ str(len(y1)) + '\n\t Percentage Actives: %.2f' % (100*sum(y1)/len(y1)))
    print ('SET 2: \n\t Num.of data points: '+ str(len(y2)) + '\n\t Percentage Actives: %.2f' % (100*sum(y2)/len(y2)))
    
    
    return x1, x2, y1, y2



def draw_covariance(data_in_columns, my_descriptors):
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    """
    Pearson's r can range from -1 to 1. An r of -1 indicates a perfect negative 
    linear relationship between variables, an r of 0 indicates no linear 
    relationship between variables, and an r of 1 indicates a perfect positive 
    linear relationship between variables.
    """

    covariance_matrix = numpy.zeros( (len(my_descriptors),len(my_descriptors)) )
    
    for k in range(len(my_descriptors)):
       
        for ki in range(len(my_descriptors)):
            a = (data_in_columns[:,k]).transpose()
            b = (data_in_columns[:,ki]).transpose()
            covariance_matrix[k, ki] = numpy.corrcoef(a, b)[0,1]
            print(my_descriptors[k], my_descriptors[ki], covariance_matrix[k, ki])
    
    
    
    fig, ax = plt.subplots()
    im = ax.imshow(covariance_matrix, cmap='YlGn', interpolation='nearest')
    #fig.colorbar(im, ax=ax)
    
    ax.set_title('Heatmap of Pearson product-moment correlation coefficient')
    print (my_descriptors)
    plt.xticks(np.arange(len(my_descriptors)), my_descriptors)
    ax.set_xticklabels(my_descriptors)
    plt.yticks(np.arange(len(my_descriptors)), my_descriptors)
    ax.set_yticklabels(my_descriptors)
    fig.colorbar(im)
    #plt.colorbar()
    plt.savefig('covar_small.png')





#----------------------------------------------------------------------------#
 
if __name__ == "__main__":
    # execute only if run as a script
    import sys
    
    # LOAD DATA
    all_x, all_y =  load_data (sys.argv[1], qm_descriptors)
    print(all_x[:,-1])
    
    # PLOT
    #draw_covariance(all_x, qm_descriptors)


    print('> DONE.')

