process Fq_QC {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0'}"
    tag { SampleName }
    publishDir "${params.Result_Folder}/${SampleName}/QC/Raw", mode: 'copy'

    input:
    tuple val(SampleName), file(reads)

    output:
    tuple val(SampleName), file("*.zip"), file("*.html")

    script:
    if (params.Seq_Tech == "Illumina") {
        """
                fastqc -memory=10000 ${reads[0]} ${reads[1]}
            """
    }
    else {
        """
                fastqc -memory=10000 ${reads}
            """
    }
}

process BAM_QC {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--hd87286a_2' :
        'biocontainers/samtools:1.17--hd87286a_2'}"
    tag { SampleName }
    publishDir "${params.Result_Folder}/${SampleName}/QC/Raw", mode: 'copy'

    input:
    tuple val(SampleName), file(alignment), file(index)

    output:
    tuple val(SampleName), file("*.tsv"), file("*.txt")

    script:
    """
        name=\$( echo ${alignment} | sed "s/.bam//g" )
        samtools flagstats -O tsv ${alignment} > \${name}_flagstats.tsv
        samtools stats ${alignment} > \${name}_stats.txt
    """
}

process Combine_QC {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'biocontainers/multiqc:1.25.1--pyhdfd78af_0'}"
    tag { SampleName }
    publishDir "${params.Result_Folder}/${SampleName}/QC", mode: 'copy'

    input:
    tuple val(SampleName), file(snpEff_csv), file(snpEff_vcf), path(raw_dat)

    output:
    tuple val(SampleName), file("${SampleName}_multiqc.html")

    script:
    """
       multiqc --filename ${SampleName}_multiqc.html Raw
    """
}

process Depths {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--hd87286a_2' :
        'biocontainers/samtools:1.17--hd87286a_2'}"
    tag { Name }
    publishDir "${params.Result_Folder}/${Name}/QC/Raw", mode: 'copy'

    input:
    tuple val(Name), file(bam), file(bai)

    output:
    tuple val(Name), file("${Name}_depths.tsv")

    script:
    """
            samtools depth -J -aa -d 1000000000 ${bam} > ${Name}_depths.tsv
        """
}

process DepthGraph {
    container 'docker://rocker/tidyverse:latest'
    tag { Name }
    publishDir "${params.Result_Folder}/${Name}/QC", mode: 'copy'

    input:
    tuple val(Name), file(depths)

    output:
    file "${Name}_depth.pdf"

    script:
    """
        PlotDepth.r ${depths} ${Name}
    """
}
