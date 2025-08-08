include {
    CombineBAM
} from '../modules/preprocess.nf'

include { fromSamplesheet } from 'plugin/nf-validation'

workflow PreProcess {
    main:
        BamFiles_ch = Channel.empty()

        Channel.fromSamplesheet("input")
            | map {row -> tuple(row.external_id[0], file(row.bam_1[0], checkIfExists: true), file(row.bam_2[0], checkIfExists: true))}
            | view()
            | set { BamFiles_ch }
        pooled_reads_ch = CombineBAM(BamFiles_ch)
    emit:
    pooled_reads = pooled_reads_ch
}
