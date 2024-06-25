

include { EXTRACT_PROMOTERS } from "../modules/local/EXTRACT_PROMOTERS.nf"
//include { SEQUENCE_SIMILARITY } from "subworkflows/SEQUENCE_SIMILARITY.nf"


workflow PROMOTER_SIMILARITY_ANALYSIS{
    
    take: 
    gtf
    fasta

    main: 

    EXTRACT_PROMOTERS(gtf, fasta)



}