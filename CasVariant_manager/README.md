# Cas Variant manager
This software automates the creation of Cas Orthologue dataset for Cloudgene. This application/pipeline does not serve any command line purpose.
Within cloudgene a dataset's access can be controled similarily to an application (restricted to certain group or individuals)
This software also helps manage the addition of new Cas Variant.

## Usage

Instalation Within cloudgene :

Step 1 : Download the CRISPR-BEasy github tool  
Step 2 : Configure the nextflow.config file (i.e. fill the cloudgene_exe and install_dir with relevant parameters)  
Step 3 : Install the cloudgene application  
```  
./cloudgene install [path]/CRISPR-BEasy/CasVariant_manager/cloudgene.yaml  
```
Step 4 : if restart cloudgene

After which the following tool can be found in runs  
![CasVariantManager](https://github.com/CERC-Genomic-Medicine/CRISPR-BEasy/blob/main/CasVariant_manager/CasVariantManager.png)

## Job Workflow
This tool creates a json file and a cloudgene.yaml (which reference the json file). The Json file contains all pertinent information to run a sgRNA_library_design and if ran via cloudgene sgRNA_library_design will automatically reference the appropriate file. See sgRNA_library_design details if ran via command-line.

1) Create a Json file and Cloudgene yaml file.  
2) Install the yaml file via the cloudgene_exe.  

