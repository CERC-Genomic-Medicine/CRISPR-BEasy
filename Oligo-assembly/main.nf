#!/usr/bin/env nextflow

params.configFile=file("./nextflow.config")
/*
========================================================================================
    CERC-Genomic-Medicine/nf-Oligomer
========================================================================================
    Github : https://github.com/CERC-Genomic-Medicine/CRISPR_Library_prep
    Author: Vincent Chapdelaine
    ---------------------------
*/

nextflow.enable.dsl = 2


/*
========================================================================================
    RUN Oligomer
========================================================================================
*/

include { Oligomer_WEB } from './module/oligomer'
include { Validate_libraries } from './module/validate_libraries'
include { Validate_Instructions } from './module/validate_instructions'
include { Validate_Target_library } from './workflows/Target.nf'
include { Validate_Negative_control_library } from './workflows/Negative_controls.nf'
include { Validate_Positive_control_library } from './workflows/Positive_controls.nf'


LP = params.Library_positive ? params.Library_positive : params.Library_positive_default
LN = params.Library_negative ? params.Library_negative : params.Library_negative_default

if(params.outdir == "default" || params.outdir == null) {
    params.pubDir = "output/${params.project}/"
} else {
    params.pubDir = "${params.outdir}/"
}

if (params.Errors == null) {
    params.Errors = params.pubDir
} else {
   params.Errors = "${params.pubDir}/Errors/"
}

if (params.Auxiliary_files == null) {
    params.Auxiliary_files = params.pubDir
} else {
   params.Auxiliary_files = "${params.pubDir}/Auxiliary_files/"
}

if (params.Oligomer_repository == null) {
    params.Oligomer_repository = params.pubDir
} else {
   params.Oligomer_repository = "${params.pubDir}/Oligomer_repository/"
}



workflow {


Channel
    .from(params.Positive_library_instruction)
    .map { content ->
        def file = File.createTempFile('Positive_instructions', '.txt')  // Create a temporary file
        file.text = content  // Write the string content to the file
        return file.path  // Return the file path to the channel 
	}
    .set { LP_instruction }  // Set the channel
    Target_ch = Channel.fromPath(params.Library_target, checkIfExists:true)
    Positive_ch = Channel.fromPath(LP, checkIfExists:true)
    Negative_ch = Channel.fromPath(LN, checkIfExists:true)

    general = Validate_libraries(Target_ch, Positive_ch, Negative_ch)
    Instruction = Validate_Instructions(LP_instruction, Positive_ch)
    if (params.Library_positive) { 
    	Positive_Annotations = Validate_Positive_control_library(Positive_ch) 
    	Annotations = Positive_Annotations.Annotation
	Oligomer_WEB(Instruction,general.CSV.collect(),Annotations)
       } else {
   		Oligomer_WEB(Instruction,general.CSV.collect(),[])
	}
}

