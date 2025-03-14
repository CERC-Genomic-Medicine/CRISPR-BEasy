

include { combine_vep } from '../modules/local/Combine_VEP'
include { combine_csv } from '../modules/local/Combine_CSV'
include { report } from '../modules/local/report'

workflow Finalization {
take :
   csv
   vep
main:
   CSV=combine_csv( csv )
   VEP=combine_vep( vep.flatten().map{ t -> [t.baseName.tokenize('_')[3], t]}.groupTuple(), CSV.remove)
   report(CSV.CSV, VEP.collect() )

}

