process SIMILARITY_MMSEQS{

    input: 
    tuple val(meta), path(fasta)
    val(promoter_length)

    output: 
    tuple val(meta), path("${meta.id}_sim.csv"), emit: sim_csv

    script:
    """
    grep -A 1 "ENSG00000162998.4" $fasta > gene.fa
    mmseqs easy-search gene.fa $fasta res tmp --prefilter-mode 1 --search-type 3 --max-seqs 1000000 --diag-score 0
    # keep columns 1 2 and 3 and add a header
    echo "gene1,gene2,promoter_similarity_mmseqs" > ${meta.id}_sim.csv
    awk -v OFS=',' '{print \$1, \$2, \$3}' res >> ${meta.id}_sim.csv
    """
}