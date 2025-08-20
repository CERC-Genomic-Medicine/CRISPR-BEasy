# sgRNA_library_design

sgRNA library design is a pipeline and tool of the CRISPR-BEasy software to facilitate the design of crispr Base Editing screens. As such the goal is to derive sgRNA from region(s) or gene(s) to create an annotated library.  

Further information about the theory behing BE sgRNA library design can be found within the [manual](https://cerc-genomic-medicine.ca/Manuals/CRISPR-BEasy/) ! 

## Installation 

Necessary software :

jq  
curl  
bowtie (version 1)  
nextflow  
apptainer
R (see [packages](https://github.com/CERC-Genomic-Medicine/CRISPR-BEasy/blob/Development/sgRNA_library_design/installed_packages_R.txt) )  
python3 (see [packages](https://github.com/CERC-Genomic-Medicine/CRISPR-BEasy/blob/Development/sgRNA_library_design/requirements_python3_env.txt) )  


### Installation For cloudgene :

Step 1 : Download the CRISPR-BEasy github tool  
Step 2 : Configure the nextflow.config file (i.e. fill the cloudgene_exe and install_dir with relevant parameters)  
Step 3 : Download the VEP ensembl sif
```  
apptainer pull --name vep.sif docker://ensemblorg/ensembl-vep
```  

Step 4 : Install the cloudgene application    
```  
./cloudgene install [path]/CRISPR-BEasy/sgRNA_library_design/cloudgene.yaml    
```

Step 5 : Restart cloudgene  

Step 6 : Populate the datasets  
1) Install the [Genome Manager](https://github.com/CERC-Genomic-Medicine/CRISPR-BEasy/tree/Development/Genome_manager) and relevant datasets  
2) Install the [Cas Variant Manager](https://github.com/CERC-Genomic-Medicine/CRISPR-BEasy/tree/Development/CasVariant_manager) and relevant datasets  

### Installation for command line

Follow step 1-3 and 6 of cloudgene installation.  

#### Regarding the datasets in command line

Within Command line datasets are refered to as Genome_Json and Cas_Variant_Json examples can be seen in assets. While within the cloudgene usage those are stored 1 genome to 1 json file, this is not necessary in command line, a concatenated file can be referenced along with the relevant id as 'genome' and 'cas_Name' parameter respectively.  
This will allow the software to extract the correct relevant information from the concatenated Genome_Json/Cas_Variant_Json. In the absence of such parameters the first entry is assumed to be correct)  
N.B. Genome_Json contains reference to other files which should be accessible on the launch (and sometimes computing ** see VEP cache **) device. AS SUCH THE DATASETS CANNOT BE TRANSFERED directly !  

## Usage Command line 
Within the command line the following parameters should be provided :  
'Genome_Json' -> json formatted file (**Files referenced should be accessible**) *  
'Cas_Variant_Json' -> json formatted file  
'Python_env'       -> Python environment containing the relevant packages  
'R_temporary_dir'  -> Where R can write temporarily (Particularily usefull in open-shared computing environment)  
'Target_genes'     -> File or \n separated list of gene or locations  
'isoform'          -> "Mane"/"None"/"Canonical"  
'GC'               -> True/False GC pattent's C should not be mutated  
'Flanking'         -> How many bp on each side of the targets should be included  
'Editors'          -> file or appropriately formatted (i.e. \t,\n) mutation instructions  
'VEP_sif'          -> ensembl's variant effect predictor's location  
  
All of these can be supplied either via --parameter value or via the nextflow.config file  








