// println "Welcome to ${params.service.name} ${workflow.manifest.version}"
// println "(c) ${workflow.manifest.author}"

include { Basic_Annotation } from '../modules/local/basic_annotation'
include { annotate } from '../Dependency/Annotation.nf'
include { to_excel } from '../modules/local/to_excel'
include { CRISPRverse } from '../modules/local/crisprVerse'
include { OnTarget } from '../modules/local/on_target'

params.Library_Type = 'Target_library'

workflow Crispr_Target_library_prep {
take :
  target_bed
  editors

main:
  target_crispr = CRISPRverse(target_bed)
  Scored = OnTarget(target_crispr.crispr_base, target_crispr.aln)
  target_BA = Basic_Annotation( Scored.scored.combine(target_bed ), editors )
  target_A = annotate( target_BA.VCF.flatten().combine(target_bed) )
  temp = target_A.Annotations
   CSV = target_BA.CSV
   output=to_excel(CSV, temp.collect(), 'Target_Library.xlsx')
emit:
   VEP = temp
   CSV = CSV

}

