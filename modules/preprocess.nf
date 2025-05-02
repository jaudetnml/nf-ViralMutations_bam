process CombineMinIONFastq {
    publishDir "${params.outdir}/${SampleName}/Combined", pattern: "${SampleName}_combined.fastq.gz", mode: 'copy'
    tag { SampleName }
    label 'process_single'

    input:
    tuple val(SampleName), path(FOLDER)

    output:
    tuple val(SampleName), path("${SampleName}_combined.fastq.gz")

    script:
    """
        if [[ "${params.Extension}" =~ .*gz ]]; then
            gunzip -c ${FOLDER}/*${params.Extension} > ${SampleName}_combined.fastq
        else
            cat ${FOLDER}/*${params.Extension} > ${SampleName}_combined.fastq
        fi
        gzip ${SampleName}_combined.fastq
    """
}

process TrimIllumina {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.24.0--heae3180_1' :
        'biocontainers/fastp:0.24.0--heae3180_1'}"
    tag { Name }
    label 'process_medium'
    publishDir "${params.outdir}/Trim", pattern: "*.fq.gz", mode: 'copy'
    publishDir "${params.outdir}/${Name}/QC/Raw", pattern: "*.json", mode: 'copy'
    publishDir "${params.outdir}/${Name}/QC", pattern: "*.html", mode: 'copy'

    input:
    tuple val(Name), file(reads1), file(reads2)

    output:
    tuple val(Name), file("${Name}_val_1.fq.gz"), file("${Name}_val_2.fq.gz"), emit: trim_reads_ch
    tuple val(Name), file("*_TrimReport.json"), file("*_TrimReport.html"), emit: trim_report_ch

    script:
    """
         fastp -i ${reads1} -I ${reads2} -o ${Name}_val_1.fq.gz -O ${Name}_val_2.fq.gz\
          -j ${Name}_TrimReport.json -h ${Name}_TrimReport.html -w ${task.cpus} ${params.TrimArgs}
    """
}

process TrimMinION {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastplong:0.2.2--heae3180_0' :
        'biocontainers/fastplong:0.2.2--heae3180_0'}"
    tag { Name }
    label 'process_medium'
    publishDir "${params.outdir}/Trim", pattern: "*.fastq.gz", mode: 'copy'
    publishDir "${params.outdir}/${Name}/QC/Raw", pattern: "*.json", mode: 'copy'
    publishDir "${params.outdir}/${Name}/QC", pattern: "*.html", mode: 'copy'

    input:
    tuple val(Name), file(reads)

    output:
    tuple val(Name), file("${Name}_trim.fastq.gz"), emit: trim_reads_ch
    tuple val(Name), file("*_TrimReport.json"), file("*_TrimReport.html"), emit: trim_report_ch

    script:
    """
        fastplong -i ${reads} -o ${Name}_trim.fastq.gz -w ${task.cpus} -m 10 -j ./${Name}_TrimReport.json -h ./${Name}_TrimReport.html  ${params.TrimArgs}
    """
}
