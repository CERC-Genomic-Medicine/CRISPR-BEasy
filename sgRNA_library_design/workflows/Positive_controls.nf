// println "Welcome to ${params.service.name} ${workflow.manifest.version}"
// println "(c) ${workflow.manifest.author}"

params.Library_Type='Positive_control_library'

include { Basic_Annotation } from '../modules/local/basic_annotation'
include { annotate } from '../Dependency/Annotation.nf'
include { to_excel } from '../modules/local/to_excel'
include { CRISPRverse } from '../modules/local/crisprVerse'
include { OnTarget } from '../modules/local/on_target'
include { split_bed } from '../modules/local/split_bed'
include { combine_general } from '../modules/local/combine_general'
include { combine_vcfs } from '../modules/local/combine_vcfs'
include { combine_annotations } from '../modules/local/combine_annotations'
include { combine_failed } from '../modules/local/combine_failed'




workflow Crispr_Positive_library_prep {
take :
  Positive_bed
  editors

main:
  Target_bed = split_bed(Positive_bed)
  //parallele
  target_crispr = CRISPRverse(Target_bed.bed_chunks.flatten())
  Scored = OnTarget(target_crispr.crispr_base, target_crispr.aln)
  target_BA = Basic_Annotation( Scored.scored.combine(Positive_bed ), editors )
  target_A = annotate( target_BA.VCF.flatten().combine(Positive_bed) )  
  //reunify
  Combine_fail = combine_failed(target_crispr.failed.collect(), 'Postitive_Control_Library')
  Combine_csv = combine_general(target_BA.CSV.collect(), 'Postitive_Control_Library') 
  vcfs = combine_vcfs(    target_BA.VCF.flatten().map( file -> [ file.getBaseName().tokenize('_')[4], file ]).groupTuple(by: [0]), 'Postitive_Control_Library',Combine_csv.ID_dic)
  vep = combine_annotations(target_A.Annotations.flatten().map( file -> [file.getBaseName().tokenize('_')[4], file]).groupTuple(by: [0]), 'Postitive_Control_Library', Combine_csv.ID_dic)
  outputCSV=Combine_csv.csv
  outputVEP=vep.tsv
  output=to_excel(Combine_csv.csv, vep.tsv.collect(),  'Postitive_Control_Library.xlsx')

emit:
  VEP = outputVEP
  CSV = outputCSV

}


