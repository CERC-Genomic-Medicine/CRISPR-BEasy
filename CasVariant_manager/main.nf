#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
import ParamDocs
import ParamUtils

/*
========================================================================================
    CERC-Genomic-Medicine/CRISPR-BEasy/sgRNA_library_design
========================================================================================
    Github : https://github.com/CERC-Genomic-Medicine/CRISPR-BEasy/tree/main/sgRNA_library_design
    Author: Vincent Chapdelaine
    ---------------------------
*/

nextflow.enable.dsl = 2


/*
========================================================================================
    RUN sgRNA_library_design
========================================================================================
*/

def paramDescriptions = ParamDocs.getParamDescriptions() // Get parameter descriptions
//=======================================================Install Database =============================================//
//------------------ Install Genome from ensembl -------------//
                def Required_parameters = ['install_cas_variant', "pam", 'protospacer_length', 'CFD_access', 'Pam_side',"rs3_access","cloudgene_exe"]
                def Ignored_parameters = ['project']
                def optional_parameters = ['install_dir','cloudgene',"cloudgene_exe"]
                def pipelineDesc = "This pipeline is a companion pipeline for CRISPR-BEasy's cloudgene data management. It creates and install a dataset corresponding to a cas variant in CRISPR-BEasy from version 1 onwards.  \n Usage nextflow main.nf --install_cas_variant <Short Name> <requirered parameters> \n  This functionnality install relevant CAS orthologue and produce a usable Json file detailling the installed CAS orthologue"
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

include { install_cas } from './modules/local/install_cas'
include { make_cloudgene } from './modules/local/cloudgene'

// ------------------------------------ MAIN WORKFLOW ----------------------------------------------------- //

workflow {
    // ================= Install Data workflow ====================== //

    if (params.install_cas_variant) {
    installed = install_cas(params.install_cas_variant)
    }
    if (params.cloudgene && params.cloudgene_exe) {
         make_cloudgene(installed.okay)
     }
}

