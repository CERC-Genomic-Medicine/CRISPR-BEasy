#!/usr/bin/env nextflow

params.configFile=file("./nextflow.config")

/*
========================================================================================
    CERC-Genomic-Medicine/nf-CRISPR_Library_prep
========================================================================================
    Github : https://github.com/CERC-Genomic-Medicine/CRISPR_Library_prep
    Author: Vincent Chapdelaine
    ---------------------------
*/

nextflow.enable.dsl = 2


/*
========================================================================================
    RUN Crispr Library Prep 
========================================================================================
*/

include { Validate } from './modules/local/validate'
include { Crispr_Target_library_prep } from './workflows/Target.nf'
include { Crispr_Negative_library_prep } from './workflows/Negative_controls.nf'
include { Crispr_Positive_library_prep } from './workflows/Positive_controls.nf'
include { Finalization } from './workflows/Combine.nf'




if(params.outdir == "default" || params.outdir == null) {
    params.pubDir = "output/${params.project}/"
} else {
    params.pubDir = "${params.outdir}/"
}


if (params.Libraries == null) {
    params.Libraries = params.pubDir
} else {
   params.Libraries = "${params.pubDir}/Libraries/"
}

if (params.Auxiliary_files == null) {
    params.Auxiliary_files = params.pubDir
} else {
   params.Auxiliary_files = "${params.pubDir}/Auxiliary_files/"
}

if (params.Report_output == null) {
    params.Report_output = params.pubDir
} else {
   params.Report_output = "${params.pubDir}/Report_output/"
}

if (params.Errors == null) {
    params.Errors = params.pubDir
} else {
   params.Errors = "${params.pubDir}/Errors/"
}




workflow {

Channel
    .from(params.Editors)
    .map { content ->
        def file = new File("${workflow.workDir}/Editors_${params.project}.txt")
        file.text = content  // Write the string content to the file
        return file.path  // Return the file path to the channel
    }
    .set { Editors_file }  // Set the channel

Channel
    .from(params.protein_id)
    .map { content ->
        def file = new File("${workflow.workDir}/Target_${params.project}.txt")
        file.text = content  // Write the string content to the file
        return file.path  // Return the file path to the channel
    }
    .set { target_ch }  // Set the channel

Channel
    .from(params.protein_Pos)
    .map { content ->
        def file = new File("${workflow.workDir}/Positive_${params.project}.txt")
        file.text = content  // Write the string content to the file
        return file.path  // Return the file path to the channel
    }
    .set { Positive_ch }  // Set the channel

Channel
    .from(params.protein_Neg)
    .map { content ->
        def file = new File("${workflow.workDir}/Negative_${params.project}.txt")
        file.text = content  // Write the string content to the file
        return file.path  // Return the file path to the channel
    }
    .set { Negative_ch }  // Set the channel

       valid = Validate(target_ch, Positive_ch, Negative_ch, Editors_file)
       Crispr_Target_library_prep( valid.target_bed, Editors_file )
       out_CSV = Crispr_Target_library_prep.out.CSV
       out_VEP = Crispr_Target_library_prep.out.VEP
       if (params.protein_Pos) { 
       		Crispr_Positive_library_prep(valid.positive_bed,Editors_file) 
       		out_CSV = out_CSV.concat(Crispr_Positive_library_prep.out.CSV)
       		out_VEP = out_VEP.concat(Crispr_Positive_library_prep.out.VEP)
       }
       if (params.protein_Neg) { 
       		Positive_ch.view()
		Crispr_Negative_library_prep(valid.negative_bed,  Editors_file) 
       		out_CSV = out_CSV.concat(Crispr_Negative_library_prep.out.CSV)
       		out_VEP = out_VEP.concat(Crispr_Negative_library_prep.out.VEP)
       }
       Finalization(out_CSV.flatten().collect(), out_VEP.flatten().collect())
}

