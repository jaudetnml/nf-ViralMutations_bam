process SetSnpEff {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.2--hdfd78af_1' :
        'biocontainers/snpeff:5.2--hdfd78af_1'}"
    label 'process_single'
    tag { snpeff_name }

    input:
    path snpeff_dataFolder
    val snpeff_name
    file snpeff_config

    output:
    file "snpEff2.config"

    script:
    """
        snpEff_line="\n${snpeff_name}.genome:${snpeff_name}"
        cat ${snpeff_config} > snpEff2.config
        
        echo \$snpEff_line >> snpEff2.config

        snpEff build -noCheckCds -noCheckProtein -c snpEff2.config -dataDir . ${snpeff_name}
    """
}