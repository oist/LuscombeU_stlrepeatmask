process SOFTMASK_STATS {
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jq:1.6':
        'biocontainers/jq:1.6' }"


    input:
    path (files)

    output:
    path "*_mqc.tsv",  emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Here we make the header
    echo "# id: 'repeat summary'" > masking_stats_mqc.tsv
    echo "# section_name: 'repeat masking summary statistics'" >> masking_stats_mqc.tsv
    echo "# format: 'tsv'" >> masking_stats_mqc.tsv
    echo "# plot_type: 'table'" >> masking_stats_mqc.tsv
    echo "# description: 'This plot shows a brief summary of each genomes whose repeats has been masked'" >> masking_stats_mqc.tsv
    echo "# pconfig:" >> masking_stats_mqc.tsv
    echo "#    id: 'repeat summary'" >> masking_stats_mqc.tsv
    echo "#    title: 'repeat summary'" >> masking_stats_mqc.tsv
    echo "#    ylab: ''" >> masking_stats_mqc.tsv
    echo "id\tTotal scaffold length\tTotal contig length\ttantan masked bases\twindowmasker masked bases\trmodeler masked bases (REPM)\trmodeler masked bases (DFAM)\trmodeler masked bases (EXTR)\tAll maskers combined" >> masking_stats_mqc.tsv
    # Here we loop on samples
    for file in *_tantan.assembly_summary
    do
        SAMPLE=\$(basename \$file _tantan.assembly_summary)
        printf "\$SAMPLE\t" >> masking_stats_mqc.tsv
        grep 'Total scaffold length' \$file | tail -n 1 | awk '{printf \$4"\t"}' >> masking_stats_mqc.tsv
        grep 'Total contig length'   \$file | tail -n 1 | awk '{printf \$4"\t"}' >> masking_stats_mqc.tsv
        grep 'soft-masked bases'     \$file | tail -n 1 | awk '{printf \$4"\t"}' >> masking_stats_mqc.tsv
        grep 'soft-masked bases' \${SAMPLE}_windowmasker.assembly_summary | awk '{printf \$4}' >> masking_stats_mqc.tsv
        for REPM_RUN in REPM DFAM EXTR
        do
            REPM_SUMMARY=\${SAMPLE}_\${REPM_RUN}.assembly_summary
            if [ -e \$REPM_SUMMARY ] ; then
                grep 'soft-masked bases' \$REPM_SUMMARY | tail -n 1 | awk '{printf "\t"\$4}' >> masking_stats_mqc.tsv
            else
                printf '\tNA' >> masking_stats_mqc.tsv
            fi
        done
        grep 'soft-masked bases' \${SAMPLE}_allmaskers.assembly_summary | awk '{printf "\t"\$4}' >> masking_stats_mqc.tsv
        printf '\n' >> masking_stats_mqc.tsv
    done
    """
}
