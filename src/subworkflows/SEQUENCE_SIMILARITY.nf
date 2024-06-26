
include { SIMILARITY_SW } from "../modules/local/SIMILARITY_SW.nf"

workflow SEQUENCE_SIMILARITY{

    take: 
    promoters_fasta
    promoter_length


    main: 

    // if(params.hamming){
    //     SIMILARITY_HAMMING(promoters_gtf, promoters_fasta)       
    // }
    matrix = Channel.fromPath("${params.matrix}")
    matrix.view()
    promoters_fasta.view()
    if(params.sw){
        promoters_fasta.view()
        SIMILARITY_SW(promoters_fasta, promoter_length, matrix)
    }


}
