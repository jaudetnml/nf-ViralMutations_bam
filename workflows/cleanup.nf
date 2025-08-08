include {
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

    main:
    depths_ch = Channel.empty()
    clean_out_ch = Channel.empty()

    Depths (aligned_reads_ch)
    | set { depths_ch }
    
    if (params.SNP_MaxCov > 0) {
        Downsample(aligned_reads_ch.join(depths_ch))
            | SI_ds
            | set { clean_out_ch }
    }
    else {
        clean_out_ch = aligned_reads_ch
    }

    emit:
    depths      = depths_ch
    final_align = clean_out_ch
}
