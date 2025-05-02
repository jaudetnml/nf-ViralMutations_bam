include {
    CombineMinIONFastq ;
    TrimMinION ;
    TrimIllumina ;
} from '../modules/preprocess.nf'

include { fromSamplesheet } from 'plugin/nf-validation'

workflow PreProcess {
    main:
    if (params.Seq_Tech == "MinION") {
        trimmed_reads_ch = Channel.empty()

        if (params.MinION_split) {
            RawReadFolders_ch = Channel.empty()

            Channel.fromSamplesheet("input")
                | map { row -> tuple(row.external_id[0], file(row.long_reads[0], type: 'dir', checkIfExists: true))}
                | view()
                | set { RawReadFolders_ch }
            CombineMinIONFastq(RawReadFolders_ch)
                | TrimMinION
                | set { trimmed_reads_ch }
        }
        else {
            MinION_reads_ch = Channel.empty()

            Channel.fromSamplesheet("input")
                | map { row -> tuple(row.external_id[0], file(row.long_reads[0], checkIfExists: true))}
                | view()
                | set { MinION_reads_ch }
            TrimMinION(MinION_reads_ch)
                | set { trimmed_reads_ch }
        }
    }
    else {
        raw_reads_ch = Channel.empty()

        Channel.fromSamplesheet("input")
                | map { row -> tuple(row.external_id[0], file(row.fastq_1[0], checkIfExists: true), file(row.fastq_2[0], checkIfExists: true))}
                | view()
                | set { raw_reads_ch }

        trimmed_reads_ch = TrimIllumina(raw_reads_ch)
    }

    emit:
    trimmed_reads = trimmed_reads_ch.trim_reads_ch
    trim_reports  = trimmed_reads_ch.trim_report_ch
}
