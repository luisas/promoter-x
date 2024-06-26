
include { SIMILARITY_SW } from "../modules/local/SIMILARITY_SW.nf"
include { SIMILARITY_MMSEQS } from "../modules/local/SIMILARITY_MMSEQS.nf"
include { ADD_PROMOTER_LENGTH } from "../modules/local/ADD_PROMOTER_LENGTH.nf"

workflow SEQUENCE_SIMILARITY{

    take: 
    promoters_fasta
    promoter_length


    main: 
    if(params.sw){
        matrix = Channel.fromPath("${params.matrix}")
        SIMILARITY_SW(promoters_fasta, promoter_length, matrix)
        sequence_similarity = SIMILARITY_SW.out.sim_csv
    }

    if(params.mmseqs){
        SIMILARITY_MMSEQS(promoters_fasta, promoter_length)
        sequence_similarity = SIMILARITY_MMSEQS.out.sim_csv
    }

    

    ADD_PROMOTER_LENGTH(sequence_similarity, promoter_length)
    csv = ADD_PROMOTER_LENGTH.out.csv

    emit: 
    csv

}
