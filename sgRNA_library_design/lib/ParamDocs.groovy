// lib/ParamDocs.groovy
class ParamDocs {
    static Map getParamDescriptions() {
        return [
    project              : "Label or identifier for the current run or project.",
    Genome_Json          : "Path to the genome reference JSON database.                                [Default Define in nextflow.config file]",

// Install parameter & Database Derived Parameter
// Tools
    Python_env           : "Location of the activate file of the Python environment                     [Default Define in nextflow.config file]",
    R_temporary          : "Location of the created temporary files / usefull for limited diskspace     [Default Define in nextflow.config file]",
    clougene_exe         : "location of the clougene executatble                                        [Default Define in nextflow.config file]",
]
}
}
