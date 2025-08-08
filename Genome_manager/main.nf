#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
import JsonUtils
import ParamDocs
import ParamUtils

/*
========================================================================================
    CERC-Genomic-Medicine/CRISPR-BEasy/Genome_Manager
========================================================================================
    Github : https://github.com/CERC-Genomic-Medicine/CRISPR-BEasy/tree/main/sgRNA_library_design
    Author: Vincent Chapdelaine
    ---------------------------
*/

nextflow.enable.dsl = 2


/*
========================================================================================
    RUN Genome_Manager
========================================================================================
*/

def paramDescriptions = ParamDocs.getParamDescriptions() // Get parameter descriptions
//=======================================================Install Database =============================================//
//------------------ Install Genome from ensembl -------------//
                def Required_parameters = ['install_ensembl']
                def Ignored_parameters = ['project']
                def optional_parameters = ['install_dir','Genome_Json','vep_threshold','cloudgene','Python_env','R_temporary_dir','cloudgene_exe']
                def pipelineDesc = "This pipeline helps design guide RNA for base editing based on genes/custom regions. \n Usage nextflow main.nf --install_ensembl <comma-delimited gneomes> \n  This functionnality install all genomes specified in the parameter from ensembl and produce a usable Json file detailling the installed genomes"
                // Validate the presence/absence of optional/Requiered/illegal parameters and produce custom errors
                ParamUtils.validateExactParamSet(
                        params,
                        Required_parameters,
                        optional_parameters,
                        ParamUtils.generateParamHelpMessage,
                        paramDescriptions,
                        pipelineDesc,
                        Ignored_parameters
                        )

// =======================================  Workflow & processess import  =================================//

include { install_workflow } from './workflow/install_workflow'
include { make_cloudgene } from './modules/local/cloudgene'

// ------------------------------------ MAIN WORKFLOW ----------------------------------------------------- //

workflow {
    // ================= Install Data workflow ====================== //

    if (params.install_ensembl) {
           Channel
            .from(params.install_ensembl.tokenize(','))
            .set { install_assemblies_ch }
	    installed=install_workflow(install_assemblies_ch)
            println('::message:: currently Dust_repeat, repeatmasquer are not implemented')
            if (params.cloudgene_exe){
                make_cloudgene(installed.cloudgene_ready)
                }
    }
}

