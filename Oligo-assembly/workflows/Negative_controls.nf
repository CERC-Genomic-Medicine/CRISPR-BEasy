// println "Welcome to ${params.service.name} ${workflow.manifest.version}"
// println "(c) ${workflow.manifest.author}"

include { Validate_Annotation } from '../module/validate_annotations'


workflow Validate_Negative_control_library {
take :
  Target_ch

main:
  Annotations = Validate_Annotation('Negative_control_library',Target_ch)
  output= Annotations.validated
emit:
  Annotation=output

}
