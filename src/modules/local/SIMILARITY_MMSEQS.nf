process SIMILARITY_MMSEQS{

    input: 
    tuple val(meta), path(fasta)
    val(column_name)

    output: 
    tuple val(meta), path("${meta.id}_sim.csv"), emit: csv

    script:
    def args = task.ext.args ?: ''
    """
    #grep -A 1 "ENSG00000162998" $fasta > gene.fa
    mmseqs easy-search $fasta $fasta res tmp $args
    echo "gene1,gene2,${column_name},aln_length,bit_score,eval" > ${meta.id}_sim.csv
    awk -v OFS=',' '{print \$1, \$2, \$3, \$4, \$5, \$11}' res >> ${meta.id}_sim.csv
    """
}