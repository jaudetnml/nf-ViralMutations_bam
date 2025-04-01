nextflow.enable.dsl = 2

include { Alignment  } from './workflows/alignment.nf'

include { PreFlight  } from './workflows/preflight.nf'

include { PreProcess } from './workflows/preprocess.nf'

include { CleanUp    } from './workflows/cleanup.nf'

include { Results    } from './workflows/results.nf'

include {
    BAM_QC ;
    DepthGraph ;
    Combine_QC
} from './modules/qc.nf'

workflow {
    setFolder = {
        def sample = it[0]
        def newPath = "${params.Result_Folder}/${sample}/QC/Raw"

        return [it[0], it[1], it[2], newPath]
    }
    
    reference_ch = Channel.fromPath(params.Target_Reference, type: 'file', checkIfExists: true)
    PreFlight()
    PreProcess()
    Alignment(PreProcess.out.trimmed_reads, PreFlight.out.Target_Reference, PreFlight.out.Host_bwa, PreFlight.out.Host_minimap, reference_ch)
    BAM_QC(Alignment.out.aligned_reads)
    CleanUp(Alignment.out.aligned_reads, PreFlight.out.Primers)
    Results(CleanUp.out.final_align, reference_ch, PreFlight.out.SnpEff_config, CleanUp.out.depths)
    DepthGraph(CleanUp.out.depths)
    Results.out.snpEff_report
        | map(setFolder)
        | Combine_QC
}
