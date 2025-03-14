// println "Welcome to ${params.service.name} ${workflow.manifest.version}"
// println "(c) ${workflow.manifest.author}"

params.Library_Type='Positive_control_library'

include { Basic_Annotation } from '../modules/local/basic_annotation'
include { annotate } from '../Dependency/Annotation.nf'
include { to_excel } from '../modules/local/to_excel'
include { CRISPRverse } from '../modules/local/crisprVerse'
include { OnTarget } from '../modules/local/on_target'

workflow Crispr_Positive_library_prep {
take :
  Positive_bed
  editors

main:
  Positive_crispr = CRISPRverse(Positive_bed)
  Scored = OnTarget(Positive_crispr.crispr_base, Positive_crispr.aln)
  Positive_BA = Basic_Annotation( Scored.scored.combine(Positive_bed ), editors )
  Positive_A = annotate( Positive_BA.VCF.flatten().combine(Positive_bed))
  temp = Positive_A.Annotations
  CSV = Positive_BA.CSV
  to_excel(CSV, temp.collect(), 'Postitive_Control_Library.xlsx')

emit:
  VEP = temp
  CSV = CSV

}


