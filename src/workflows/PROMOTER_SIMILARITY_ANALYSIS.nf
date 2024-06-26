

include { EXTRACT_PROMOTERS } from "../modules/local/EXTRACT_PROMOTERS.nf"
include { SEQUENCE_SIMILARITY } from "../subworkflows/SEQUENCE_SIMILARITY.nf"


workflow PROMOTER_SIMILARITY_ANALYSIS{
    
    take: 
    input_files
    promoter_length

    main: 

    EXTRACT_PROMOTERS(input_files, promoter_length)

    SEQUENCE_SIMILARITY(EXTRACT_PROMOTERS.out.promoters_fa, promoter_length)

     

}