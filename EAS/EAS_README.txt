Requirements
	rdkit/2018.03.1-foss-2017a-Python-3.6.6
	Python/3.6.6-foss-2017a
	gaussian/16.a.03-avx2
	numpy/1.13.3-foss-2017a-Python-3.6.6
	scikit-learn/0.19.2-foss-2017a-Python-3.6.6
	gepol93/8
	chargemol/3.5
Paths & global variables
	All paths and global variables are stored in the file “MAP_paths.py”.
	This files also contains templates that are used to generate input 
	files and submission scripts (to run jobson cluster).
	Depending on how the call to the cluster will be made, some parts of 
	this file may be removed/changed.

To run
ex:
	python path_to_MAP/MAP.py compute some_sdf.sdf
	python path_to_MAP/MAP.py compute "COc1cccnc1"
This will run the calculations to generate all the QM descriptors required. 
You may want to create a dedicated folder to run in, since the program uses 
some generic file names. Once the gaussian jobs are done, you can run the 
following to extract the descriptors and run the random forest model:
	python path_to_MAP.py extract 


