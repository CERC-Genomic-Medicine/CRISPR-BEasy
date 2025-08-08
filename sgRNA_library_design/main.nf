#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
import JsonUtils
import ParamDocs
import ParamUtils

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

def paramDescriptions = ParamDocs.getParamDescriptions() // Get parameter descriptions
// ==================================== Documentation & Loading params ==========================//
if (CLOUDGENE_WORKSPACE_TYPE ) {
//------------------- Run Cloudgene ------------//
                def Required_parameters = ['Genome_Json','Cas_Variant_Json', 'Python_env', 'R_temporary_dir', 'Target_genes', 'isoform', 'GC', 'Flanking', 'Editors', 'feature', 'VEP_sif']
                def Ignored_parameters = ['send_mail',"project","remove_empty","output","outdir","Libraries","Auxiliary_files","Report_output","Report_ancillary","Errors","Library_Type","configFile","config-file","pubDir","pub-dir"]
                def optional_parameters = ['R_temporary_dir', 'Positive_genes', 'Negative_genes','CFD_Threshold', 'CFD_Count','limit_bp','soft_bp_limit','CFD_files','crispr_chunksize']
                def pipelineDesc = "This pipeline helps design guide RNA for base editing based on genes/custom regions. \n Usage of this version should be restricted to cloudgene_backend"
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
                def asmParams = JsonUtils.loadAssemblyParams(params.Genome_Json, null)
                params.putAll(asmParams)
                def casParams = JsonUtils.loadCasVariant(params.Cas_Variant_Json, null)
                params.putAll(casParams)
} else {
//---------------------------- Run locally -----------//
                def Required_parameters = ['Genome_Json','Cas_Variant_Json', 'Python_env', 'R_temporary_dir', 'Target_genes', 'isoform', 'GC', 'Flanking', 'Editors', 'feature', 'VEP_sif']
                def Ignored_parameters = ['send_mail','genome','cas_Name',"output","outdir","Libraries","Auxiliary_files","Report_output","Report_ancillary","Errors","Library_Type","configFile","config-file","pubDir","pub-dir"]
                def optional_parameters = ['Python_env','R_temporary_dir','Positive_genes', 'Negative_genes','CFD_Threshold', 'CFD_Count','limit_bp','soft_bp_limit','CFD_files','crispr_chunksize']
                def pipelineDesc = "This pipeline helps design guide RNA for base editing based on genes/custom regions. \n Usage : nextflow main.nf --Genome <value> [--Genome_Json <value> ] --Editors <value> --cas_Name <value> [--Cas_Variant_Json <value>] --Target_genes <value> [--Positive_genes <value>] [--Negative_genes <value>] [--feature_type <value>] [--Isoform <value>] [--GC <value>] [--Flanking <value>] [--CFD_Theshold <value> --CFD_count <value>] [--CFD_files <value>] [--VEP_sif <value>] \n"
                ParamUtils.validateExactParamSet(
                // Validate the presence/absence of optional/Requiered/illegal parameters and produce custom errors
                        params,
                        Required_parameters,
                        optional_parameters,
                        ParamUtils.generateParamHelpMessage,
                        paramDescriptions,
                        pipelineDesc,
                        Ignored_parameters
                        )
                def asmParams = JsonUtils.loadAssemblyParams(params.Genome_Json, params.genome)
                asmParams.each { k, v -> params[k] = v }
                def casParams = ParamUtils.loadCasVariant(params.CasVariantJson, params.cas_Name)
                casParams.each { k, v -> params[k] = v }

}



// ----------------------- Load Process and workflow -------------------------//

include { Validate } from './modules/local/validate'
include { Crispr_Target_library_prep } from './workflows/Target.nf'
include { Crispr_Target_library_prep as Crispr_Negative_library_prep } from './workflows/Target.nf'
include { Crispr_Target_library_prep as Crispr_Positive_library_prep } from './workflows/Target.nf'
include { Finalization } from './workflows/Combine.nf'


//============================================ Main Workflow ==================================//
workflow {

//------------------------- Load Files ---------------------//

// Targets //
Target_name = Channel.value('Study_Target')
Channel
    .from(params.Target_genes)
    .map { content ->
        def inputFile = new File(content.toString())
        if (inputFile.exists()) {
            return inputFile.path  // Use the existing file
        } else {
            def file = new File("${workflow.workDir}/Target_${params.project}.txt")
            file.text = content.toString()  // Write string content to new file
            return file.path  // Return the new file path
        }
    }
    .set { Target_ch }

// Postive Controls //
Positive_name = Channel.value('Positive_controls')
Channel
    .from(params.Positive_genes ?: "")
    .map { content ->
        def path = content.toString()
        def inputFile = new File(path)

        if (inputFile.exists()) {
            return inputFile.path  // Use the existing file
        }

        def file = new File("${workflow.workDir}/Positive_${params.project}.txt")
        file.text = path ? path : ""  // Write string or leave empty
        return file.path
    }
    .set { Positive_ch }


// Negative Controls //
Negative_name = Channel.value('Negative_controls')
Channel
    .from(params.Negative_genes ?: "")
    .map { content ->
        def path = content.toString()
        def inputFile = new File(path)

        if (inputFile.exists()) {
            return inputFile.path  // Use the existing file
        }

        def file = new File("${workflow.workDir}/Negative_${params.project}.txt")
        file.text = path ? path : ""  // Write string or leave empty
        return file.path
    }
    .set { Negative_ch }

// Editors //
Channel
    .from(params.Editors)
    .map { content -> 
        def inputFile = new File(content.toString())
        if (inputFile.exists()) {  
            return inputFile.path  // Use the existing file
        } else {
            def file = new File("${workflow.workDir}/Editor_${params.project}.txt")
            file.text = content.toString()  // Write string content to new file
            return file.path  // Return the new file path
        }
    }
    .set { Editors_file }

// GFF database (derived from Genome) //
Channel
    .fromPath(params.gffdb_loc).set { GFF_database }

// Fasta (derived from Genome) //
Channel
    .fromPath(params.fasta_loc)
    .map { path -> file(path.toString() + '.fai') }
    .set { fasta_database }

// Boyle_lab Blacklist //
Channel
    .from(params.Boyle_Lab ?: "")
    .map { content ->
        def path = content.toString()
        def inputFile = new File(path)

        if (inputFile.exists()) {
            return inputFile.path  // Use the existing file
        }

        def file = new File("${workflow.workDir}/Boyle_lab_${params.project}.bed.gz")
        file.text = ""  // Write string or leave empty
        return file.path
    }
    .set { boyle_ch }

// ------------------------------ Process ----------------------------------//

       valid = Validate(Target_ch, Positive_ch, Negative_ch, Editors_file, GFF_database, fasta_database, boyle_ch)
        Crispr_Target_library_prep( valid.target_bed, Editors_file, Target_name )
       out_CSV = Crispr_Target_library_prep.out.CSV
       out_VEP = Crispr_Target_library_prep.out.VEP

        if (params.Positive_genes) { 
       		Crispr_Positive_library_prep(valid.positive_bed, Editors_file, Positive_name )
       		out_CSV = out_CSV.concat(Crispr_Positive_library_prep.out.CSV)
       		out_VEP = out_VEP.concat(Crispr_Positive_library_prep.out.VEP)
       }
       if (params.Negative_genes) { 
		Crispr_Negative_library_prep(valid.negative_bed, Editors_file, Negative_name ) 
       		out_CSV = out_CSV.concat(Crispr_Negative_library_prep.out.CSV)
       		out_VEP = out_VEP.concat(Crispr_Negative_library_prep.out.VEP)
       }
       Finalization(out_CSV.flatten().collect(), out_VEP.flatten().collect())

}

workflow.onComplete {
    def user_email = System.getenv('CLOUDGENE_USER_EMAIL')
    def SMTP_USER = System.getenv('CLOUDGENE_SMTP_USER')
    def service_url = System.getenv('CLOUDGENE_SERVICE_URL')
    def service_name = System.getenv('CLOUDGENE_SERVICE_NAME')
    def user_full_name = System.getenv('CLOUDGENE_USER_FULL_NAME')
    def SMTP_PASS = System.getenv('CLOUDGENE_SMTP_PASSWORD')
    System.setProperty("mail.smtp.auth", "true")
    System.setProperty("mail.smtp.starttls.enable", "true")
    System.setProperty("mail.smtp.host", "smtp-mail.outlook.com")
    System.setProperty("mail.smtp.port", "587")
    System.setProperty("mail.smtp.user", SMTP_USER)
    System.setProperty("mail.smtp.password", SMTP_PASS)

    println(CLOUDGENE_SERVICE_URL)
    //job failed
    if (!workflow.success) {
        def statusMessage = workflow.exitStatus != null ? "failed" : "canceled"
        if ( params.send_mail && user_email != null){
            sendMail{
                to "${user_email}"
                subject "[${service_name}] Job ${params.project} failed."
                body "Dear ${user_full_name}, your ${service_name} job has ${statusMessage}. More details can be found at the following link: https://crispr-beasy.cerc-genomic-medicine.ca/#!jobs/${params.project}"

            }
        }
        println "::error:: ${service_name} Server job ${statusMessage}."
        return
    }  else {
        //job successful
                println " Attempt ::message:: Sent email notification to <b> $user_email </b> from ${SMTP_USER}"
                sendMail{
                to "${user_email}"
                subject "[${service_name}] Job ${params.project} is complete"
                body "Dear ${user_full_name}, \n your ${service_name} job has finished succesfully.\n You can download the results from the following link: https://crispr-beasy.cerc-genomic-medicine.ca/index.html#!jobs/${params.project}"
                }
}
}
