process CombineMinIONFastq {
    publishDir "${params.Result_Folder}/${SampleName}/Combined", pattern: "${SampleName}_combined.fastq.gz", mode: 'copy'
    tag { SampleName }

    input:
    tuple val(SampleName), val(Barcode)

    output:
    tuple val(SampleName), path("${SampleName}_combined.fastq.gz")

    script:
    """
        if [[ "${params.Extension}" =~ .*gz ]]; then
            gunzip -c ${params.Data_Folder}/${Barcode}/*${params.Extension} > ${SampleName}_combined.fastq
        else
            cat ${params.Data_Folder}/${Barcode}/*${params.Extension} > ${SampleName}_combined.fastq
        fi
        gzip ${SampleName}_combined.fastq
    """
}

process TrimIllumina {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trim-galore:0.6.10--hdfd78af_1' :
        'biocontainers/trim-galore:0.6.10--hdfd78af_1'}"
    tag { Name }
    label 'TrimIllumina'
    publishDir "${params.Result_Folder}/Trim", pattern: "*.fq.gz", mode: 'copy'
    publishDir "${params.Result_Folder}/${Name}/QC/Raw", pattern: "*.txt", mode: 'copy'

    input:
    tuple val(Name), file(reads)

    output:
    tuple val(Name), file("${Name}_val_1.fq.gz"), file("${Name}_val_2.fq.gz"), emit: trim_reads_ch
    tuple val(Name), file("*_trimming_report.txt"), emit: trim_report_ch

    script:
    """
        trim_galore -q 25 --gzip --max_n 5 --paired -j 4 ${params.Illumina_TrimArgs}--basename ${Name} ${reads[0]} ${reads[1]}
    """
}

process TrimIlluminaClip {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trim-galore:0.6.10--hdfd78af_1' :
        'biocontainers/trim-galore:0.6.10--hdfd78af_1'}"
    tag { Name }
    label 'TrimIllumina'
    publishDir "${params.ResultFolder}/Trim", pattern: "*.fq.gz", mode: 'copy'
    publishDir "${params.ResultFolder}/${Name}/QC/Raw", pattern: "*.txt", mode: 'copy'

    input:
    tuple val(Name), file(reads)

    output:
    tuple val(Name), file("${Name}_val_1.fq.gz"), file("${Name}_val_2.fq.gz"), emit: trim_reads_ch
    tuple val(Name), file("*_trimming_report.txt"), emit: trim_report_ch

    script:
    """
        trim_galore -q 25 --gzip --max_n 5 --paired -j 4 --clip_R1 ${params.Illumina_clipR1} --clip_R2 ${params.Illumina_clipR2} --basename ${Name} ${reads[0]} ${reads[1]}
    """
}

process TrimIlluminaCustAdapt {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trim-galore:0.6.10--hdfd78af_1' :
        'biocontainers/trim-galore:0.6.10--hdfd78af_1'}"
    tag { Name }
    label 'TrimIllumina'
    publishDir "${params.ResultFolder}/Trim", pattern: "*.fq.gz", mode: 'copy'
    publishDir "${params.ResultFolder}/${Name}/QC/Raw", pattern: "*.txt", mode: 'copy'

    input:
    tuple val(Name), file(reads)

    output:
    tuple val(Name), file("${Name}_val_1.fq.gz"), file("${Name}_val_2.fq.gz"), emit: trim_reads_ch
    tuple val(Name), file("*_trimming_report.txt"), emit: trim_report_ch

    script:
    """
        trim_galore -q 25 --gzip --max_n 5 -a ${params.Illumina_adapt1} -a2 ${params.Illumina_adapt2} --paired -j 4 --basename ${Name} ${reads[0]} ${reads[1]}
    """
}

process TrimMinION {
    cpus '16'
    memory '8G'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/porechop:0.2.3_seqan2.1.1--0' :
        'biocontainers/porechop:0.2.3_seqan2.1.1--0'}"
    tag { Name }
    publishDir "${params.Result_Folder}/Trim", pattern: "*.fastq.gz", mode: 'copy'
    publishDir "${params.Result_Folder}/${Name}/QC/Raw", pattern: "*.txt", mode: 'copy'

    input:
    tuple val(Name), file(reads)

    output:
    tuple val(Name), file("${Name}_trim.fastq.gz"), emit: trim_reads_ch
    tuple val(Name), file("*porechop_log.txt"), emit: trim_report_ch

    script:
    """
        porechop -i ${reads} -o ${Name}_trim.fastq.gz -t ${task.cpus} --check_reads 1000 > ${Name}_porechop_log.txt
    """
}
