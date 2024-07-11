process CUSTOMMODULE {
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jq:1.6':
        'biocontainers/jq:1.6' }"


    input:
    path (assemt)
    path (assemw)
    path (assemr)

    output:
    path "*_mqc.tsv",  emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    echo "# id: 'repeat summary'" > repeatmqc_mqc.tsv
    echo "# section_name: 'repeat masking summary statistics'" >> repeatmqc_mqc.tsv
    echo "# format: 'tsv'" >> repeatmqc_mqc.tsv
    echo "# plot_type: 'bargraph'" >> repeatmqc_mqc.tsv
    echo "# description: 'This plot shows a brief summary of each genomes whose repeats has been masked'" >> repeatmqc_mqc.tsv
    echo "# pconfig:" >> repeatmqc_mqc.tsv
    echo "#    id: 'repeat summary'" >> repeatmqc_mqc.tsv
    echo "#    title: 'repeat summary'" >> repeatmqc_mqc.tsv
    echo "#    ylab: ''" >> repeatmqc_mqc.tsv
    echo "id\tTotal scaffold length\tTotal contig length\ttantan masked bases\twindowmasker masked bases\trmodeler masked bases" >> repeatmqc_mqc.tsv
    printf "\$(basename $assemt .assembly_summary)\t" >> repeatmqc_mqc.tsv
    grep 'Total scaffold length' $assemt | tail -n 1 | awk '{print \$4}' | tr '\n' '\t' >> repeatmqc_mqc.tsv
    grep 'Total contig length' $assemt | tail -n 1 | awk '{print \$4}' | tr '\n' '\t' >> repeatmqc_mqc.tsv
    grep 'soft-masked bases' $assemt | tail -n 1 | awk '{print \$4}' | tr '\n' '\t' >> repeatmqc_mqc.tsv
    grep 'soft-masked bases' $assemw | tail -n 1 | awk '{print \$4}' | tr '\n' '\t' >> repeatmqc_mqc.tsv
    grep 'soft-masked bases' $assemr | tail -n 1 | awk '{print \$4}' >> repeatmqc_mqc.tsv


    echo "# id: 'tantan repeat summary'" > tantanmqc_mqc.tsv
    echo "# section_name: 'tantan repeat masking summary statistics'" >> tantanmqc_mqc.tsv
    echo "# format: 'tsv'" >> tantanmqc_mqc.tsv
    echo "# plot_type: 'heatmap'" >> tantanmqc_mqc.tsv
    echo "# description: 'This plot shows a brief summary of each genomes whose repeats has been masked'" >> tantanmqc_mqc.tsv
    echo "# pconfig:" >> tantanmqc_mqc.tsv
    echo "#    id: 'tantan repeat summary'" >> tantanmqc_mqc.tsv
    echo "#    title: 'tantan repeat summary'" >> tantanmqc_mqc.tsv
    echo "#    ylab: ''" >> tantanmqc_mqc.tsv
    echo "id\tTotal scaffold length\tTotal contig length\tsoft masked bases" >> tantanmqc_mqc.tsv
    for i in $assemt
    do
        printf "\$(basename \$i .assembly_summary)\t" >>tantanmqc_mqc.tsv
        grep 'Total scaffold length' \$i | tail -n 1 | awk '{print \$4}' | tr '\n' '\t' >> tantanmqc_mqc.tsv
        grep 'Total contig length' \$i | tail -n 1 | awk '{print \$4}' | tr '\n' '\t' >> tantanmqc_mqc.tsv
        grep 'soft-masked bases' \$i | tail -n 1 | awk '{print \$4}' >> tantanmqc_mqc.tsv
    done


    echo "# id: 'windowmasker repeat summary'" > windowmqc_mqc.tsv
    echo "# section_name: 'windowmasker repeat masking summary statistics'" >> windowmqc_mqc.tsv
    echo "# format: 'tsv'" >> windowmqc_mqc.tsv
    echo "# plot_type: 'heatmap'" >> windowmqc_mqc.tsv
    echo "# description: 'This plot shows a brief summary of each genomes whose repeats has been masked'" >> windowmqc_mqc.tsv
    echo "# pconfig:" >> windowmqc_mqc.tsv
    echo "#    id: 'windowmasker repeat summary'" >> windowmqc_mqc.tsv
    echo "#    title: 'windowmasker repeat summary'" >> windowmqc_mqc.tsv
    echo "#    ylab: ''" >> windowmqc_mqc.tsv
    echo "id\tTotal scaffold length\tTotal contig length\tsoft masked bases" >> windowmqc_mqc.tsv
    for i in $assemw
    do
        printf "\$(basename \$i .assembly_summary)\t" >> windowmqc_mqc.tsv
        grep 'Total scaffold length' \$i | tail -n 1 | awk '{print \$4}' | tr '\n' '\t' >> windowmqc_mqc.tsv
        grep 'Total contig length' \$i | tail -n 1 | awk '{print \$4}' | tr '\n' '\t' >> windowmqc_mqc.tsv
        grep 'soft-masked bases' \$i | tail -n 1 | awk '{print \$4}' >> windowmqc_mqc.tsv
    done


    echo "# id: 'repeatmasker repeat summary'" > repeatmaskmqc_mqc.tsv
    echo "# section_name: 'repeatmasker repeat masking summary statistics'" >> repeatmaskmqc_mqc.tsv
    echo "# format: 'tsv'" >> repeatmaskmqc_mqc.tsv
    echo "# plot_type: 'heatmap'" >> repeatmaskmqc_mqc.tsv
    echo "# description: 'This plot shows a brief summary of each genomes whose repeats has been masked'" >> repeatmaskmqc_mqc.tsv
    echo "# pconfig:" >> repeatmaskmqc_mqc.tsv
    echo "#    id: 'repeatmasker repeat summary'" >> repeatmaskmqc_mqc.tsv
    echo "#    title: 'repeatmasker repeat summary'" >> repeatmaskmqc_mqc.tsv
    echo "#    ylab: ''" >> repeatmaskmqc_mqc.tsv
    echo "id\tTotal scaffold length\tTotal contig length\tsoft masked bases" >> repeatmaskmqc_mqc.tsv
    for i in $assemr
    do
        printf "\$(basename \$i .assembly_summary)\t" >> repeatmaskmqc_mqc.tsv
        grep 'Total scaffold length' \$i | tail -n 1 | awk '{print \$4}' | tr '\n' '\t' >> repeatmaskmqc_mqc.tsv
        grep 'Total contig length' \$i | tail -n 1 | awk '{print \$4}' | tr '\n' '\t' >> repeatmaskmqc_mqc.tsv
        grep 'soft-masked bases' \$i | tail -n 1 | awk '{print \$4}' >> repeatmaskmqc_mqc.tsv
    done
    """
}
