process Sort_Index {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--hd87286a_2' :
        'biocontainers/samtools:1.17--hd87286a_2'}"
    tag { Name }
    label 'process_low'
    publishDir "${params.outdir}/${Name}/Alignments", mode: 'copy'

    input:
    tuple val(Name), file(samfile)

    output:
    tuple val(Name), file("*.bam"), file("*.bam.bai")

    script:
    """
        filename=\$( echo ${samfile} | cut -f 1 -d '.' )
        samtools view -u ${samfile} | samtools sort -@ ${task.cpus} -o \${filename}_sorted.bam
        samtools index \${filename}_sorted.bam
    """
}
