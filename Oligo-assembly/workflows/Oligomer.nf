// println "Welcome to ${params.service.name} ${workflow.manifest.version}"
// println "(c) ${workflow.manifest.author}"
if(params.outdir == "default" || params.outdir == null) {
    params.pubDir = "output/${params.project}/Oligomer/"
} else {
    params.pubDir = "${params.outdir}/Oligomer/"
}

include { oligomer } from '../modules/local/oligmer'


workflow Oligomer {
take:
  Library_ch
  Annotation_ch
  

main:
  oligomer(Library_ch,Annotation_ch)
}


// workflow.onComplete {
//
//    def report = new CloudgeneReport()

    //job failed
//    if (!workflow.success) {
//        if (params.config.send_mail){
 //           sendMail{
//                to "${params.user.email}"
//                from "${params.service.email}"
//                subject "[${params.service.name}] Job ${params.project} failed."
//                body "Hi ${params.user.name}, your job failed :(. Logs can be accessed at ${params.service.url}/index.html#!jobs/${params.project}"
//
//            }
//        }
//        report.error("Job failed. Reason: " + workflow.errorMessage)
//        return
//    } else {
//        report.ok("Job terminated successfully. Duration: " + workflow.duration)
//    }
//
    //job successful
//    if (params.config.send_mail){
//        sendMail{
//            to "${params.user.email}"
//            from "${params.service.email}"
//            subject "[${params.service.name}] Job ${params.project} is complete."
//            body "Hi ${params.user.name}, your job completed successfully and can be accessed at ${params.service.url}/index.html#!jobs/${params.project}"
 //       }
//        report.ok("Sent email notification to <b>${params.user.email}</b>")
//    } else {
//
//    }
//}
