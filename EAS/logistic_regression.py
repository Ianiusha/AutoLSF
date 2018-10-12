# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 14:27:51 2018


Author: Anna Tomberg

This code is used to train and test a logistic regression model 
with chosen parameters.
Example of usage:
    
    python logistic_regression.py train my_data.csv
    
"""

# LOAD LIBRARIES
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler 
from sklearn.model_selection import GridSearchCV
from sklearn.feature_selection import RFE

from sklearn.metrics import accuracy_score, roc_auc_score, f1_score
from sklearn.metrics import classification_report,confusion_matrix
import data_preprocess, numpy
from data_preprocess import qm_descriptors




def scan_parameters(x,y):
 
    x = numpy.delete(x, len(qm_descriptors), axis=1)

    penalties = ['l1', 'l2']
    Cs = [0.001, 0.01, 0.1, 1, 10]
    weights = [None, 'balanced']
    
    params_dict = {'C': Cs, 'class_weight' : weights, 'penalty': penalties}
    
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    scaler.fit(x) 
    x = scaler.transform(x) 
    
    grid_search = GridSearchCV(LogisticRegression(solver='liblinear'), params_dict, cv = 3, n_jobs = 12)
    grid_search.fit(x, y)
    print(grid_search.best_score_ , grid_search.best_params_)
    
    return grid_search.best_params_


def perform_RFE(x,y):
    
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    scaler.fit(x) 
    x = scaler.transform(x) 
    
    estimator = LogisticRegression(solver='liblinear')
    selector = RFE(estimator, 1, step=1)
    selector = selector.fit(x, y)
    
    return selector.ranking_


def train_logist(X_train, y_train, X_valid, y_valid, param_list = None):
    
    x_train = numpy.delete(X_train, len(qm_descriptors), axis=1)
    x_valid = numpy.delete(X_valid, len(qm_descriptors), axis=1)
    
    
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    scaler.fit(x_train)
    
    x_train = scaler.transform(x_train)  
    
    if param_list == None:
        logreg = LogisticRegression(solver='liblinear')
        
    else:
        logreg = LogisticRegression(solver='liblinear', C = param_list['C'], penalty = param_list['penalty'], class_weight = param_list['class_weight'])

        
    
    logreg.fit(x_train, y_train)
    
    
    x_valid = scaler.transform(x_valid)  
    predictions = logreg.predict(x_valid)
    #probabilities = clf.predict_proba(x_valid) 
    
        
    #print(confusion_matrix(y_valid,predictions))
    print(classification_report(y_valid,predictions)) 
    accuracy = (100*accuracy_score(y_valid, predictions))    
    print('Accuracy: %.2f' % accuracy)
    
    
    return accuracy





def build_logist(filename, option, model_name = None):  
     
    
    # LOAD DATA
    descriptors = qm_descriptors
    X, Y = data_preprocess.load_data(filename, descriptors)

    
    if option == 'default':
        
        print('Training Logist...')
        print('*-----------------------------*')
        print('Training on default parameters.')
        
        accuracies_default = []
        for i in range(10):
            x_train, x_valid, y_train, y_valid = data_preprocess.split_data (X, Y, partition=0.20)
            accuracies_default.append(train_logist(x_train, y_train, x_valid, y_valid))
        
        print('Average accuracy over 10 default runs: %.2f' % numpy.mean(accuracies_default))
        
        
    elif option == 'train':
        
        print('*-----------------------------*')
        print('Searchig for best parameters.')
        
        params = []
        accuracies = []

        
        for i in range(10):
            x_train, x_valid, y_train, y_valid = data_preprocess.split_data (X, Y, partition=0.20)
            best_parameters = scan_parameters(x_train, y_train)
            params.append(best_parameters)
            accuracy = train_logist(x_train, y_train, x_valid, y_valid, best_parameters)
            accuracies.append(accuracy)
            
        print('*-----------------------------*')
        print('Summary of Results.')            
        print('*-----------------------------*')
        
        
        for i in range(len(accuracies)):
            print('Run ' + str (i+1)+ ' ', params[i], ' : ', accuracies[i])
        
    elif option == 'RFE':
        
        print('*-----------------------------*')
        print('Recursive feature estimation.')
        #http://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.RFE.html#sklearn.feature_selection.RFE
        
        ranking = perform_RFE(X, Y)
            
        print('*-----------------------------*')
        print('Ranking of descriptors.')            
        print('*-----------------------------*')
        for d in range(len(qm_descriptors)):
            print(qm_descriptors[d], ranking[d])
        
  
        

    elif option == 'test':
        
        print('TESTING')
        print('*-----------------------------*')
        #penalties = 'l2'
        #Cs = 0.001
        #weights = None
        
        penalties = 'l1'
        Cs = 10
        weights = None
        
        
        params_dict = {'C': Cs, 'class_weight' : weights, 'penalty': penalties}
        print(params_dict)
        
        acc_list = []
        for i in range(10):
            x_train, x_valid, y_train, y_valid = data_preprocess.split_data (X, Y, partition=0.20)
            acc_list.append(train_logist(x_train, y_train, x_valid, y_valid,params_dict))
        
        print('Summary of Results.')            
        print('*-----------------------------*')
        print('Average accuracy over 20 runs: %.2f' % numpy.mean(acc_list))
        
        


        
if __name__ == "__main__":
    # execute only if run as a script
    import sys
    
    

    option = sys.argv[1]
    filename = sys.argv[2]
    
    if len(sys.argv) > 3:
        model_n = sys.argv[3]
    else:
        model_n = None
    
    build_logist(filename, option, model_n)

    print('> DONE.')

