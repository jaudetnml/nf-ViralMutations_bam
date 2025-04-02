include {
    CombineMinIONFastq ;
    TrimMinION ;
    TrimIllumina ;
} from '../modules/preprocess.nf'

workflow PreProcess {
    main:
    if (params.Seq_Tech == "MinION") {
        if (params.MinION_split) {
            channel.fromPath("${params.MinION_Samples}", type: 'file')
                | splitCsv(header: true)
                | map { row -> tuple(row.Sample_Name, row.Barcode)}
                | view()
                | set { RawReadFolders_ch }
            CombineMinIONFastq(RawReadFolders_ch)
                | TrimMinION
                | set { trimmed_reads_ch }
        }
        else {
            channel.fromPath("${params.Data_Folder}/*${params.Extension}", type: 'file')
                | set { MinION_reads_ch }
            TrimMinION(MinION_reads_ch)
                | set { trimmed_reads_ch }
        }
    }
    else {
        channel.fromFilePairs(["${params.Data_Folder}/**_R{1,2}_001${params.Extension}", "${params.Data_Folder}/**_R{1,2}${params.Extension}"], maxDepth: 2, syntax: "glob", type: "file")
            | set { raw_reads_ch }

        trimmed_reads_ch = TrimIllumina(raw_reads_ch)
    }

    emit:
    trimmed_reads = trimmed_reads_ch.trim_reads_ch
    trim_reports  = trimmed_reads_ch.trim_report_ch
}
