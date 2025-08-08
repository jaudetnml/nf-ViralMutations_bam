include {
    SetSnpEff;
} from '../modules/preflight.nf'

workflow PreFlight{
    main:
        snpeff_folder_ch = channel.fromPath("${params.SnpEff_Folder}/${params.SnpEff_Name}", type: 'dir')
        snpeff_config_ch = channel.fromPath("${projectDir}/data/snpEff.config")
        SetSnpEff(snpeff_folder_ch, params.SnpEff_Name, snpeff_config_ch)

    emit:
        SnpEff_config = SetSnpEff.out.ifEmpty("EMPTY")