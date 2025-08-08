// lib/ParamDocs.groovy
class ParamDocs {
    static Map getParamDescriptions() {
        return [
    project              : "Label or identifier for the current run or project.",
    install_dir          : "Directory where resources and data will be installed.                      [Default Define in nextflow.config file]",
    cloudgene            : "Whether to install the dataset as an app in a cloudgene instance [Flag] [default: False]",

// Install parameter & Database Derived Parameter
    install_cas_variant : 'The name of the cas orthologue to install',
    pam                 : 'The sequence of the protospacer adjacent motif',
    protospacer_length  : 'The length of the protospacer',
    CFD_access          : 'Whether it should be able to access to CDF scoring, warning message for all but SpCas9 is hard coded [True/False]',
    Pam_side            : 'Which side the pam should be found                                           [3prime/5prime]',
    rs3_access          : 'Whether it should be able to access rs3 scoring, warning message for all but SpCas9 is hard coded [True/False]',
// Tools
    clougene_exe         : "location of the clougene executatble                                        [Default Define in nextflow.config file]"
]
}

