// println "Welcome to ${params.service.name} ${workflow.manifest.version}"
// println "(c) ${workflow.manifest.author}"

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

params.Library_Type = 'Target_library'

workflow Crispr_Target_library_prep {
take :
  target_bed
  editors

main:
    //parallele
  Target_bed = split_bed(target_bed)
  target_crispr = CRISPRverse(Target_bed.bed_chunks.flatten())
  Scored = OnTarget(target_crispr.crispr_base, target_crispr.aln)
  target_BA = Basic_Annotation( Scored.scored.combine(target_bed ), editors )
  target_A = annotate( target_BA.VCF.flatten().combine(target_bed) )  
  //Reunify
  Combine_csv = combine_general(target_BA.CSV.collect(), 'Study_Target_library') 
Combine_fail = combine_failed(target_crispr.failed.collect(),  'Study_Target_library')
  vcfs = combine_vcfs(    target_BA.VCF.flatten().map( file -> [ file.getBaseName().tokenize('_')[4], file ]).groupTuple(by: [0]), 'Study_Target_library', Combine_csv.ID_dic)
  vep = combine_annotations(target_A.Annotations.flatten().map( file -> [file.getBaseName().tokenize('_')[4], file]).groupTuple(by: [0]), 'Study_Target_library', Combine_csv.ID_dic)
  output=to_excel(Combine_csv.csv, vep.tsv.collect(), 'Target_Library.xlsx')
  outputCSV=Combine_csv.csv
  outputVEP=vep.tsv
emit:
   VEP = outputVEP
   CSV = outputCSV

}

