process CombineBAM{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--hd87286a_2' :
        'biocontainers/samtools:1.17--hd87286a_2'}"
    tag { Name }
    label 'process_single'

    publishDir "${params.outdir}/${Name}/Alignments", mode: 'copy'

    input:
    tuple val(Name), path(bam1), path(bam2)

    output:
    tuple val(Name), path("${Name}_combined.bam"), path("${Name}_combined.bam.bai")

    script:
    """
            samtools merge -o ${Name}_combined.bam ${bam1} ${bam2}
            samtools index ${Name}_combined.bam
    """
}