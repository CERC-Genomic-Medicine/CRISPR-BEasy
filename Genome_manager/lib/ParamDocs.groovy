// lib/ParamDocs.groovy
class ParamDocs {
    static Map getParamDescriptions() {
        return [
    install_ensembl      : "Install a genome from Ensembl by assembly name.\n    Usage: --install_ensembl <genome>",
    project              : "Label or identifier for the current run or project.",
    install_dir          : "Directory where resources and data will be installed.                      [Default Define in nextflow.config file]",
    vep_threshold        : "VEP cache size threshold (in GiB) to decide whether to symlink or extract the tar.gz.     [Default Define in nextflow.config file]",
    cloudgene            : "Whether to install the dataset as an app in a cloudgene instance [Flag] [default: False]",

// Tools
    Python_env           : "Location of the activate file of the Python environment                     [Default Define in nextflow.config file]",
    R_temporary          : "Location of the created temporary files / usefull for limited diskspace     [Default Define in nextflow.config file]",
    clougene_exe         : "location of the clougene executatble                                        [Default Define in nextflow.config file]",
]
}
}
