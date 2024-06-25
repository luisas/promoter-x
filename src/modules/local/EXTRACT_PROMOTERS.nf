process EXTRACT_PROMOTERS{

    input: 
    tuple val(meta), path(gtf), path(fasta)
    val(promoter_length)

    output:
    tuple val(meta), path("promoters_${promoter_length}.gtf"), emit: promoters_gtf
    tuple val(meta), path("promoters_${promoter_length}.fa"), emit: promoters_fa


    script: 
    """
    01_extract_promoters.py ${gtf} ${fasta} promoters_${promoter_length}.gtf promoters_${promoter_length}.fa ${promoter_length}
    """
}