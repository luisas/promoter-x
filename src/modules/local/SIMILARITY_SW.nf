process SIMILARITY_SW{


    input: 
    tuple val(meta), path(fasta)
    val(promoter_length)
    path(matrix)

    output: 
    tuple val(meta), path("*_sim.csv"), emit: sim_csv

    script:
    """
    head -n 10 $fasta > TEST.fasta
    swps3 $matrix TEST.fasta TEST.fasta -j 1 > sws.csv
    sed -i 's/x/\t/g' sws.csv
    sw_to_sim.py sws.csv $promoter_length ${meta.id}_sim.csv
    """
}