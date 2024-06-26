process SIMILARITY_SW{

    input: 
    tuple val(meta), path(fasta)
    val(promoter_length)
    path(matrix)

    output: 
    tuple val(meta), path("${meta.id}_sim.csv"), emit: sim_csv

    script:
    """
    # select line with gene id = gene id and corresponding fasta
    grep -A 1 "ENSG00000162998.4" $fasta > gene.fa
    swps3 $matrix gene.fa $fasta -j 1 > sws.csv
    sed -i 's/x/\t/g' sws.csv
    sw_to_sim.py sws.csv $promoter_length ${meta.id}_sim.csv
    """
}