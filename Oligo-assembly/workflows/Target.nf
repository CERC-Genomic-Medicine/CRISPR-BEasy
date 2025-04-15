// println "Welcome to ${params.service.name} ${workflow.manifest.version}"
// println "(c) ${workflow.manifest.author}"

include { Validate_Annotation } from '../module/validate_annotations'


workflow Validate_Target_library {
take :
  Target_ch

main:
  Annotations = Validate_Annotation('Target_library',Target_ch)
  output= Annotations.validated
emit:
  Annotation=output

}
