include {
    Dedup ;
    Remove_secondaries ;
    PrimerClip ;
    Downsample
} from '../modules/cleanup.nf'

include {
    Sort_Index as SI_ds
} from '../modules/alignment.nf'

include {
    Depths
} from '../modules/qc.nf'

workflow CleanUp {
    take:
    aligned_reads_ch
    primers_ch

    main:
    depths_ch = Channel.empty()
    clean_out_ch = Channel.empty()
    if (params.Primer_Locs) {
        if (params.Seq_Tech == "Illumina") {
            aligned_reads_ch
                | Dedup
            PrimerClip(Dedup.out.dedup_data_ch.combine(primers_ch))
                | Remove_secondaries
                | Depths
                | set { depths_ch }
        }
        else {
            aligned_reads_ch.combine(primers_ch)
                | PrimerClip
                | Remove_secondaries
                | Depths
                | set { depths_ch }
        }
    } else {
        if (params.Seq_Tech == "Illumina") {
            Dedup(aligned_reads_ch)
            Remove_secondaries(Dedup.out.dedup_data_ch)
                | Depths
                | set { depths_ch }
        }
        else {
            Remove_secondaries(aligned_reads_ch)
                | Depths
                | set { depths_ch }
        }
    }
    if (params.SNP_MaxCov > 0) {
        Downsample(Remove_secondaries.out.join(depths_ch))
            | SI_ds
            | set { clean_out_ch }
    }
    else {
        clean_out_ch = Remove_secondaries.out
    }

    emit:
    depths      = depths_ch
    final_align = clean_out_ch
}
