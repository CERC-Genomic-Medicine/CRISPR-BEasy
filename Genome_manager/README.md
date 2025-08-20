# Genome Variant manager
This software automates the creation of genome data compatible with CRISPR-BEasy via the ensembl ftp repository.  
Within cloudgene a dataset's access can be controled similarily to an application (i.e. restricted to certain group or individuals).  
For command-line purpose, this pipeline (main.nf/nextflow.config) allow the creation of the files necessary for sgRNA_library_design.

## Usage

### Requiered software

jq  
curl  
bowtie (version 1)  
nextflow  
R (see [packages](https://github.com/CERC-Genomic-Medicine/CRISPR-BEasy/blob/main/sgRNA_library_design/installed_packages_R.txt))  
python3 (see [packages](https://github.com/CERC-Genomic-Medicine/CRISPR-BEasy/blob/main/sgRNA_library_design/requirements_python3_env.txt) )  

### Instalation For cloudgene :

Step 1 : Download the CRISPR-BEasy github tool  
Step 2 : Configure the nextflow.config file (i.e. fill the cloudgene_exe and install_dir with relevant parameters)  
Step 3 : Install the cloudgene application  
```  
./cloudgene install [path]/CRISPR-BEasy/Genome_manager/cloudgene.yaml  
```

Step 4 : if restart cloudgene

After which the following tool can be found in runs  
![GenomeManager](https://github.com/CERC-Genomic-Medicine/CRISPR-BEasy/blob/main/Genome_manager/GenomeManager.png)

Step 5 (optional but highly recommanded) : Restrict usage to administrator only.

### Cloudgene Usage :

Genome names are only compatible with ensembl (e.g. GRCh38.p14) see [ensembl website](https://useast.ensembl.org/index.html)

### Usage for Comandline

Follow step 1 and 2 of cloudgene installation.

Genomes can be installed via the command-line :
```  
nextflow run main.nf -c nextflow.config --install_ensembl [Desired Ensembl ID]
```

### Outputs 

The main outputs will be a json file with relevant data and data location and a cloudgene yaml file (can be ignored if running in command line)  
** if  running in command line ** the json files can be concatenated (in a json safe manner) to be referenced later as Genome_Json.

Fasta, Bsgenome, gff, gffutils database and VEP cache will also be downloaded in the install_dir location.

## Workflow

The workflow is separated in 3 process 

1) Download of relevant file *Requieres Internet*
      * includes the creation of a gffutils database
2) Creation of a bowtie index
3) Creation of a R BSgenome package
4) (Depending on option) For large VEP cache, unpacking of the tar.gz file
      * Depending on the computing environment copying and decompressing a large tar.gz file might be unfeasable.
