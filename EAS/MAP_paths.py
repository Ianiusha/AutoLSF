# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 10:57:52 2018

@author: kntb373

* In this file, the paths to helper programs are collected. 
* The keyword for each program recognizable by MAP is followed
* by the path to it (as a string in quotation marks). You must 
* make sure that all the paths are correct before launching MAP.py.

* OTHER CONSTANTS USED THROUGHOUT THE PROGRAM ARE ALSO SET HERE.

* ALL THE TEMPLATES ARE SET HERE AS WELL. (EX: template for gaussian input)


"""
from string import Template



" PATHS TO HELPER PROGRAMS "
GAUSSIAN="aussian/16.b.01-avx2"

GEPOL93="gepol93/8" # executable  = sas.out

CHARGEMOL_PATH="chargemol_09_26_2017"

RANDOM_FOREST_MODEL = "my_RF_model3_eq.sav"


" DEFAULT VALUES FOR VARIABLES "

ALL_PROPERTIES = ['Q', 'CH_Bond_order', 'Bond_orders_sum', 'SAS' , 'Fukui','CH_rSQM_E']


"""
* Van der Waals radii library
* from http://periodictable.com/Properties/A/VanDerWaalsRadius.v.html
* published and maintained by Mathematica (Wolfram)
"""
VDW_DICTIONARY = {'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'P': 1.8, 'S': 1.9, 'Cl':1.75, 'Br': 1.85, 'I': 1.98, 'B': 2.0, 'Si': 2.1}
DEFAULT_VDW_RADIUS = 2.0		# IN SAS CALCULATION (if an atom is not defines in VDW_DICTIONARY)

CONVERT_AU2KCAL_MOL = 627.5
ENERGY_CUTOFF = 1.0
ENERGY_CUTOFF2 = 3.0

NODE_MEMORY = 648               # total memory on cluster node that will be divided by the number of cpus
DEFAULT_number_of_CPUs = 2
DEFAULT_TIME="04:00:00"

MAX_CONF=20					# IN SMILES2GAUSSIAN
MINI_Iterations = 200			# IN SMILES2GAUSSIAN



SAS_INP = Template("""TITL= Calculation of the Accessible surface
TITL= $NAME
ASURF
RSOL=$SOLVENT_RAD
NDIV=5
COOF=$NAME.XYZR
SPHF=$NAME.SPH

""" )

SAS_XYZR = Template("""* Calculation of the Accessible surface                                          
* Coordinates for $NAME                
      $ATOM_NUMBER
""" )



# USE THIS TO MAKE A JOB CONTROL FILE FOR  DDEC6 CHARGES CALCULATION
DDEC_JOB_CONTROL = Template("""<atomic densities directory complete path>
$PATH_TO_CHARGEMOL/atomic_densities/
</atomic densities directory complete path>

<input filename>
$TITLE.wfx
</input filename>

<charge type>
DDEC6
</charge type>

<compute BOs>
.true.
</compute BOs>

""" )

# USE THIS TO MAKE SUBMIT FILE FOR  DDEC6 CHARGES CALCULATION
DDEC_JOB_SUBMIT = Template(\
"""#!/bin/bash
#SBATCH --workdir=$PATH
#SBATCH --nodes=1
#SBATCH --hint=nomultithread
#SBATCH --time=04:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=1000
#SBATCH --job-name=$NAME

source /etc/bashrc
$PATH_TO_CHARGEMOL/chargemol_FORTRAN_09_26_2017/compiled_binaries/linux/Chargemol_09_26_2017_linux_parallel job_control.txt
                            
""" )

# USE THIS TO RUN AN OPTIMIZATION 
OPT = Template("""%nprocshared=12
%mem=1500MB
%Chk=$TITLE.chk
#n  OPT TPSSh/Def2SVP

 $TITLE

$CHARGE $MULT
$XYZ

--Link1--
%Chk=$TITLE.chk
%NoSave
#p TPSSh/Def2TZVP GFINPUT IOP(6/7=3) 6d density=current output=wfx Geom=Check Guess=Read
    
    $TITLE - wfx

0 1

$TITLE.wfx


""" )

# USE THIS TO PREPARE OUTPUT FOR REGIOSQM (PM3)
PM3 = Template("""%nprocshared=4
%mem=1500MB
%Chk=$TITLE.chk
#p pm3 opt=(maxcycles=200) scrf=(smd,solvent=chloroform) 

 $TITLE

$CHARGE $MULT
$XYZ
--Link1--
%Chk=$TITLE.chk
%NoSave
#p pm3 Geom=Check Guess=Read scrf=(check,smd,solvent=chloroform,read,externaliteration,dovacuum) 
    
TITLE.com

$CHARGE $MULT


""" )


GAUSSIAN_SUBMIT = Template("""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --hint=nomultithread
#SBATCH --time=$TIME
#SBATCH --ntasks=$CPU
#SBATCH --mem-per-cpu=500
#SBATCH --job-name=$NAME

source /etc/bashrc
module load gaussian/16.b.01-avx2

gaussian/16.b.01-avx2/g16/g16 < ./$NAME.com > ./$LOG.log
# Loaded GAS PACKAGES:
""")


##############################################################################


ALL_OPTIONS = { 'OPT': OPT, 'DDEC_JOB_CONTROL':DDEC_JOB_CONTROL,\
               'DDEC_JOB_SUBMIT': DDEC_JOB_SUBMIT, 'PM3': PM3}

##############################################################################
