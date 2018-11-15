# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 10:29:53 2018

@author: kntb373
"""


#import numpy as np
import re,os, random, subprocess
import MAP_paths


def parse_SDF(sdf_file_name):
    
    try:
        fo = open(sdf_file_name, 'r')
        data = fo.readlines()
        fo.close()
    except:
        print (">> ERROR: cannot open sdf file.")
        return ""
    
    xyz_block = ''
    
    # find and crop XYZ block
    for i in range (len(data)-1):
        m = re.search('((-)?(\d)+.(\d)+(\s){1,}){3}([a-zA-Z]{1,2})', data[i])
        if m:
            # we found the first line in the coordinates block
            # need to save the charge and multiplicity
            t = (m.group(0)).split()
            coord = t[-1] + '\t' + t[0]+ '\t' + t[1]+ '\t' + t[2]+ '\n'
            xyz_block = xyz_block + (coord)
            #print (coord)

    return xyz_block


def parse_gaussianLOG(log_file_name):
    
    try:
        fo = open(log_file_name, 'r')
        data = fo.readlines()
        fo.close()
    except:
        print (">> ERROR: cannot open sdf file.")
        return 
    
    
    xyz_block = ''
    charge = -1000
    out_block = ""
    optimized = False
    # find and crop XYZ block
    for (i, line) in enumerate(data):
       
       if re.search("1\\\\1\\\\", line):
           optimized = True
           for line in (data[i:]):   
               out_block = out_block + (line.rstrip()).replace("\n", "")

           break   
    if optimized :
        out_block = out_block.replace(" ", "")
        # the xyz block will the in the 4th position:
        temp = out_block.split('\\\\')
        temp = (temp[3]).split('\\')
        #print(temp)
    
        charge = (temp[0].split(','))[0]
        multiplicity = (temp[0].split(','))[1]
        
        #xyz_block = xyz_block + ("%s, %s \n") % (charge,multiplicity)
        
        for i in (temp[1:]) :
            #print (i)
            atom = (i.split(','))[0]
            x = (i.split(','))[-3]
            y = (i.split(','))[-2]
            z = (i.split(','))[-1]
            
            xyz_block = xyz_block + ("%s \t %s \t %s \t %s \n") % (atom, x,y,z)
            #print "%s: %s" % (key, value)
    
        #print (xyz_block)
        return (charge,multiplicity, xyz_block)

    else :
        print (">> ERROR: log file not complete.")
        return 


def create_gaussianINP(filename, coordinates_list, charge, multiplicity, option):
    
    my_template = MAP_paths.ALL_OPTIONS.get(option)
    my_input = (my_template.substitute(TITLE=(filename.split('.com')[0]), CHARGE=charge, MULT=multiplicity, XYZ=coordinates_list))
    
    fw = open(filename, 'w')
    fw.write(my_input)
    fw.close()
    
    return (filename)



def create_DDEC_control(filename):
    
    
    name1 = 'job_control.txt'
    n = filename.split('.com')[0]
    my_template = MAP_paths.ALL_OPTIONS.get('DDEC_JOB_CONTROL')
    my_input = (my_template.substitute(TITLE=n , PATH_TO_CHARGEMOL = MAP_paths.CHARGEMOL_PATH ))
    
    fw = open(name1, 'w')
    fw.write(my_input)
    fw.close()
    
    name2 = 'submit_DDEC.sh'  
    my_template = MAP_paths.ALL_OPTIONS.get('DDEC_JOB_SUBMIT')
    my_input = (my_template.substitute(PATH = os.getcwd(), NAME = filename.split('.com')[0], PATH_TO_CHARGEMOL = MAP_paths.CHARGEMOL_PATH) )
    
    fw = open(name2, 'w')
    fw.write(my_input)
    fw.close()
    #print('>> ', filename)
    
    return name1,name2



def submit_g09_job(g09_input, num_CPUs=None, runtime = None, logNAME = None):
    
     job = 'NAN' 
    
     com_NAME = g09_input.split('.com')[0]
     
     if num_CPUs is None:
         num_CPUs = MAP_paths.DEFAULT_number_of_CPUs
         
     if runtime is None:
         runtime = MAP_paths.DEFAULT_TIME
         
     if logNAME is None:
         logNAME = com_NAME    
    
     sub = 'submit'+str(random.uniform(1, 100000))
            
     fw = open(sub, 'w')
     fw.write(MAP_paths.GAUSSIAN_SUBMIT.substitute(NAME = com_NAME, CPU = num_CPUs, TIME = runtime, LOG = logNAME))
     fw.close()
            
     
     try:
         job =  subprocess.check_output(['sbatch', sub]) 
         #job=sub
         print (job)
     except subprocess.CalledProcessError as exc:
         print(">> Error in running Gaussian job:", exc.returncode, exc.output)  


     return job
    


def MAP_generate_g09(my_sdf, charge, option):
    
    """
        Here we will ...
            - extract the xyz coordinates of the sdf file given
            - produce g09 input files for the jobs listed in options
            
    """    
    # 1. get the xyz coordinate in sdf file
    coordinates = parse_SDF(my_sdf)    
    
    # 2. create the correct name for gaussian input file
    gaussian_name = my_sdf.split('.sdf')[0]+'_'+option+'.com'   
        
    # 3. prepare and submit gaussian input file(s) to compute data for descriptors

    file_created = create_gaussianINP(gaussian_name, coordinates, charge, 1, option)
    print('>> ', file_created)
    
    if option == 'OPT':
        CPUs = 12
    else:
        CPUs = None
        
    sub_out = submit_g09_job(file_created, num_CPUs=CPUs)  
    job_id = (str(sub_out).split()[-1])[:-3]
    #job_id = 'NAN'
    if not(job_id == 'NAN') and option == "OPT":
        
    # 4. for some options, other files are also required, might as well make them now
        control, subscript = create_DDEC_control(file_created)
        try:
            dependency = ("--dependency=afterok:"+job_id)
            job =  subprocess.check_output(['sbatch', dependency, subscript]) 
            #print (job)

        except subprocess.CalledProcessError as exc:
            print(">> Error in running Gaussian job:", exc.returncode, exc.output)  

       
        

    return file_created
    
    





 
