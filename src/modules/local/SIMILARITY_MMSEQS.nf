process SIMILARITY_MMSEQS{

    input: 
    tuple val(meta), path(fasta)
    val(column_name)

    output: 
    tuple val(meta), path("${meta.id}_sim.csv"), emit: csv

    script:
    def args = task.ext.args ?: ''
    """
    grep -A 1 "ENSG00000162998" $fasta > gene.fa
    mmseqs easy-search gene.fa $fasta res tmp --max-seqs 1000000 $args
    echo "gene1,gene2,${column_name}" > ${meta.id}_sim.csv
    awk -v OFS=',' '{print \$1, \$2, \$3}' res >> ${meta.id}_sim.csv
    """
}