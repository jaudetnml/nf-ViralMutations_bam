include {
    CombineMinIONFastq ;
    TrimMinION ;
    TrimIllumina ;
} from '../modules/preprocess.nf'

workflow PreProcess {
    main:
    if (params.Seq_Tech == "MinION") {
        trimmed_reads_ch = Channel.empty()

        if (params.MinION_split) {
            RawReadFolders_ch = Channel.empty()

            channel.fromPath("${params.input}", type: 'file')
                | splitCsv(header: true)
                | map { row -> tuple(row.sample, file(row.fastq_1, type: 'dir', checkIfExists: true))}
                | view()
                | set { RawReadFolders_ch }
            CombineMinIONFastq(RawReadFolders_ch)
                | TrimMinION
                | set { trimmed_reads_ch }
        }
        else {
            MinION_reads_ch = Channel.empty()

            channel.fromPath("${params.input}", type: 'file')
                | splitCsv(header: true)
                | map { row -> tuple(row.sample, file(row.fastq_1, checkIfExists: true))}
                | view()
                | set { MinION_reads_ch }
            TrimMinION(MinION_reads_ch)
                | set { trimmed_reads_ch }
        }
    }
    else {
        raw_reads_ch = Channel.empty()

        channel.fromPath("${params.input}", type: 'file')
                | splitCsv(header: true)
                | map { row -> tuple(row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true))}
                | view()
                | set { raw_reads_ch }

        trimmed_reads_ch = TrimIllumina(raw_reads_ch)
    }

    emit:
    trimmed_reads = trimmed_reads_ch.trim_reads_ch
    trim_reports  = trimmed_reads_ch.trim_report_ch
}
