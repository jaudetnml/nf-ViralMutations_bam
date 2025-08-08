process Downsample {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-bb5f3dab55f89ca6e9acdff899d8409efffcc444:949930df72decfcefd14c5d64ef58250f319589e-0' :
        'biocontainers/mulled-v2-bb5f3dab55f89ca6e9acdff899d8409efffcc444:949930df72decfcefd14c5d64ef58250f319589e-0'}"
    tag { Name }
    label 'process_single'
    publishDir "${params.outdir}/${Name}/Alignments", mode: 'copy'

    input:
    tuple val(Name), file(bam), file(bai), file(depth)

    output:
    tuple val(Name), file("${Name}_Aligned_ds.bam")

    script:
    if (params.Seq_Tech == "Illumina") {
        """
            DownsampleToCoverage.py -r ${bam} -n ${Name}_Aligned_ds.bam -c ${depth} -t ${params.SNP_MaxCov}
        """
    }
    else {
        """
            DownsampleToCoverage.py -r ${bam} -n ${Name}_Aligned_ds.bam -c ${depth} -t ${params.SNP_MaxCov}
        """
    }
}


