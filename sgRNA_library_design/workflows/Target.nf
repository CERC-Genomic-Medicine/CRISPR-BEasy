// println "Welcome to ${params.service.name} ${workflow.manifest.version}"
// println "(c) ${workflow.manifest.author}"

include { Basic_Annotation } from '../modules/local/basic_annotation'
include { annotate } from '../modules/local/Annotation'
include { annotate_tar } from '../modules/local/Annotation_tar'
include { to_excel } from '../modules/local/to_excel'
include { CRISPRverse } from '../modules/local/crisprVerse'
include { OnTarget } from '../modules/local/on_target'
include { split_bed } from '../modules/local/split_bed'
include { combine_general } from '../modules/local/combine_general'
include { combine_vcfs } from '../modules/local/combine_vcfs'
include { combine_annotations } from '../modules/local/combine_annotations'
include { combine_failed } from '../modules/local/combine_failed'

params.Library_Type = 'Target_library'

workflow Crispr_Target_library_prep {
take :
  target_bed
  editors
  name_output
  ensembl_filter

main:

if (params.vep_cache_loc.endsWith('.tar.gz')) {
        Channel
             .fromPath(params.vep_cache_loc)
             .set { Vep_cache_ch }
}
// bsgenome_loc

Channel
    .fromPath(params.bsgenome_loc)
    .set { BSgenome_ch }


def BowtiedirPath = new File(params.bowtie_index_loc).getParent()
Channel
    .fromPath("${BowtiedirPath}", type: 'dir')
    .set { bowtie_index_folder_ch }
// Fasta (derived from Genome) //
Channel
    .fromPath(params.fasta_loc)
    .map { path -> tuple (path, file(path.toString() + '.fai')) }
    .set { fasta_database_ch }

Channel
    .fromPath("${params.CFD_files}", type: 'dir')
    .set { CFD_file_ch }
  //parallele
  Target_bed = split_bed(target_bed)
  target_crispr = CRISPRverse(Target_bed.bed_chunks.flatten(),  BSgenome_ch, CFD_file_ch, bowtie_index_folder_ch)
  Scored = OnTarget(target_crispr.crispr_base, target_crispr.aln)
  
  input_target = Scored.scored
    .combine(target_bed)
  target_BA = Basic_Annotation( input_target.combine(fasta_database_ch), editors )
  vcfs = target_BA.VCF.flatten().combine(target_bed)
  if (!params.vep_cache_loc.endsWith('.tar.gz')) {
        target_A = annotate(vcfs.combine(fasta_database_ch),ensembl_filter)
  } else {
        target_A = annotate_tar(vcfs.combine(fasta_database_ch),
    Vep_cache_ch,ensembl_filter
)

  }
  //Reunify

  Combine_csv = combine_general(target_BA.CSV.collect(), name_output.map { it + '_library' }) 
Combine_fail = combine_failed(target_crispr.failed.collect(),  name_output.map { it + '_library' })
  vcfs = combine_vcfs(    target_BA.VCF.flatten().map( file -> [ file.getBaseName().tokenize('_')[4], file ]).groupTuple(by: [0]), name_output.map { it + '_library' } , Combine_csv.ID_dic)
  vep = combine_annotations(target_A.Annotations.flatten().map( file -> [file.getBaseName().tokenize('_')[4], file]).groupTuple(by: [0]), name_output.map { it + '_library' }, Combine_csv.ID_dic)
  output=to_excel(Combine_csv.csv, vep.tsv.collect(), name_output.map { it + '_library' })

  outputCSV=Combine_csv.csv
  outputVEP=vep.tsv
emit:
   VEP = outputVEP
   CSV = outputCSV

}

