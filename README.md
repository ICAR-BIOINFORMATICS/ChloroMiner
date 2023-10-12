# ChloroMiner (Chloroplast genome extraction and assembly Tool)
The automation script to assemble chloroplast genome from whole genome NGS (forward and reverse reads).
“ChloroMiner” the chloroplast genome assembly Pipeline is a user-friendly tool designed to facilitate the assembly of chloroplast genomes from the raw sequencing data. This manual provides a step-by-step guide on how to use the pipeline effectively. 


### Before starting please ensure that you have the following dependencies installed before proceeding:


Python (version 3 or higher):https://docs.python-guide.org/starting/install3/linux/


Trimmomatic: https://git.launchpad.net/ubuntu/+source/trimmomatic/log/?h=ubuntu/focal

## Download the ChloroMiner package from the github using 

		https://github.com/ICAR-BIOINFORMATICS/ChloroMiner.git 


 
 
## After downloading the pipeline package, you need to follow the steps below one by one. 
 
Step 1: Downloaded python scripts to a desired location on your computer 

		automation_script.py and c_miner.py

Step 2: Preparing the Input Data Before running the pipeline, make sure you have the following input data ready:
		
		Raw sequencing reads in farstq.gz format.


Running the Pipeline Follow these steps to run the ChloroMiner:

1: Open a terminal or command prompt on your computer (Ctrl+Alt+T).

2: Navigate to the directory where the pipeline is installed.

3: Execute the pipeline script using the following command:

		python3 automatioon_script.py -adap <Trimomatic_adapter_path> -ref <reference_chloroplast_genome_file> -threads <no of threads>

Here -adap is the path of trimmomatic adapters in your computer; -ref is the refrence Chloroplast genome; -threads is no of threads as per your computer's specification i.e. 8/12/16/31/64


## Developed by:

Dr. Samarth Godara, Scientist, Division of Computer Applications, ICAR-Indian Agricultural Statistics Research Institute (IASRI), Library Avenue, Pusa, New Delhi-110012 For any query, contact: samarth.godara@gmail.com.

Dr. Shbana Begam, Scientist, ICAR-National Institute for Plant Biotechnology, New Delhi-110012



