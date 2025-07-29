// lib/ParamDocs.groovy
class ParamDocs {
    static Map getParamDescriptions() {
        return [
    install_ensembl      : "Install a genome from Ensembl by assembly name.\n    Usage: --install_ensembl <genome>",
    update_genomes       : "Update all genomes in the database to the most recent Ensembl release.\n    Usage: --update_genomes",
    install_url          : "Install a genome manually from provided URLs.\n    Usage: --install_url <genome> --fasta_url <url> --gff3_url <url> [--vep_cache_url <url>]",
    install_file         : "Install a genome manually from provided files \n   Usage: --install_file <genome> fasta_loc <location> --gff3_loc <location> [--vep_cache_url <url>]",
    make                 : "Install all genomes listed in the referenced genome database.\n    Usage: --make",

    project              : "Label or identifier for the current run or project.",
    install_dir          : "Directory where resources and data will be installed.                      [Default Define in nextflow.config file]",
    vep_threshold        : "VEP cache size threshold (in GiB) to decide whether to symlink or extract the tar.gz.     [Default Define in nextflow.config file]",
    Genome_Json          : "Path to the genome reference JSON database.                                [Default Define in nextflow.config file]",
    cloudgene            : "Whether to install the dataset as an app in a cloudgene instance [Flag] [default: False]",

// Install parameter & Database Derived Parameter
    release              : "Genome release version.",
    source               : "Annotation source (e.g., Ensembl, RefSeq).",
    mane_select          : "MANE annotation status (e.g., MANE_Select).",
    name                 : "Organism name (used for naming the BSgenome).",

    vep_cache_url        : "URL for the VEP cache tar.gz (indexed format).",
    fasta_url            : "URL for the genome FASTA file (.fa.gz).",
    gff3_url             : "URL for the genome GFF3 annotation file (.gz).",
    Boyle_Lab_url        : "URL for Boyle Lab blacklist regions (.bed.gz).",
    Dust_repeat_url      : "URL for Dust repeat regions (.bed.gz).",
    repeatmasquer_url    : "URL for RepeatMasker regions (.bed.gz).",

    fasta_loc            : "Local path to the genome FASTA file and its index (.fai).",
    gffdb_loc            : "Local path to the gffutils database.",
    gff3_loc             : "Local path to the GFF3 annotation file.",
    vep_cache_loc        : "Local path to the VEP cache tar.gz.",
    bowtie_index_loc     : "Prefix path to the Bowtie index files.",
    bsgenome_loc         : "Path to the BSgenome R package folder.",
    Boyle_Lab_loc        : "Path to local Boyle Lab blacklist file.",
    Dust_repeat_loc      : "Path to local Dust repeat file.",
    repeatmasquer_loc    : "Path to local RepeatMasker file.",
// Tools
    Python_env           : "Location of the activate file of the Python environment                     [Default Define in nextflow.config file]",
    R_temporary          : "Location of the created temporary files / usefull for limited diskspace     [Default Define in nextflow.config file]",
    clougene_exe         : "location of the clougene executatble                                        [Default Define in nextflow.config file]",
]
}
}
