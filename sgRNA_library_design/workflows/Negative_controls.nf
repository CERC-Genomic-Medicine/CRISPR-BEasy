// println "Welcome to ${params.service.name} ${workflow.manifest.version}"
// println "(c) ${workflow.manifest.author}"

params.Library_Type='Negative_control_library'

include { Basic_Annotation } from '../modules/local/basic_annotation'
include { annotate } from '../Dependency/Annotation.nf'
include { to_excel } from '../modules/local/to_excel'
include { CRISPRverse } from '../modules/local/crisprVerse'
include { OnTarget } from '../modules/local/on_target'

workflow Crispr_Negative_library_prep {
take :
  Negative_bed
  editors

main:
  Negative_crispr = CRISPRverse(Negative_bed)
  Scored = OnTarget(Negative_crispr.crispr_base, Negative_crispr.aln)
  Negative_BA = Basic_Annotation( Scored.scored.combine(Negative_bed ), editors )  
  Negative_A = annotate( Negative_BA.VCF.flatten().combine(Negative_bed) )
  temp = Negative_A.Annotations
  CSV = Negative_BA.CSV
  to_excel(CSV, temp.collect(), 'Negative_Control_Library.xlsx')

emit:
  VEP = temp
  CSV = CSV


}


