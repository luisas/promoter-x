

include { EXTRACT_PROMOTERS } from "../modules/local/EXTRACT_PROMOTERS.nf"
include { SEQUENCE_SIMILARITY } from "../subworkflows/SEQUENCE_SIMILARITY.nf"
include { ADD_EXPRESSION      } from "../modules/local/ADD_EXPRESSION.nf"


workflow PROMOTER_SIMILARITY_ANALYSIS{
    
    take: 
    input_files
    promoter_length
    expressions

    main: 

    EXTRACT_PROMOTERS(input_files, promoter_length)
    SEQUENCE_SIMILARITY(EXTRACT_PROMOTERS.out.promoters_fa, promoter_length)
    ADD_EXPRESSION(SEQUENCE_SIMILARITY.out.csv, expressions)


}