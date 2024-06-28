

include { EXTRACT_PROMOTERS } from "../modules/local/EXTRACT_PROMOTERS.nf"
include { SEQUENCE_SIMILARITY } from "../subworkflows/SEQUENCE_SIMILARITY.nf"
include { ADD_EXPRESSION      } from "../modules/local/ADD_EXPRESSION.nf"
include { EXTRACT_PROTEIN_SEQUENCES  } from "../modules/local/EXTRACT_PROTEIN_SEQUENCES.nf"
include { SIMILARITY_MMSEQS as SIMILARITY_MMSEQS_PROTEIN } from "../modules/local/SIMILARITY_MMSEQS.nf"
include { MERGE_DATAFRAMES } from "../modules/local/MERGE_DATAFRAMES.nf"

workflow PROMOTER_SIMILARITY_ANALYSIS{
    
    take: 
    input_files
    promoter_length
    expressions

    main: 

    translated_fasta = Channel.fromPath("${params.translated_fasta}")
    // Extract the similarity of protein sequences
    EXTRACT_PROTEIN_SEQUENCES(input_files, translated_fasta)
    SIMILARITY_MMSEQS_PROTEIN(EXTRACT_PROTEIN_SEQUENCES.out.fasta, "similarity_protein")
    
    // // Extract promoters and compare them
    // EXTRACT_PROMOTERS(input_files, promoter_length)
    // SEQUENCE_SIMILARITY(EXTRACT_PROMOTERS.out.promoters_fa, promoter_length)

    // MERGE_DATAFRAMES(SEQUENCE_SIMILARITY.out.csv, SIMILARITY_MMSEQS_PROTEIN.out.csv)

    // Add expression data
    //ADD_EXPRESSION(MERGE_DATAFRAMES.out.csv, expressions)
    ADD_EXPRESSION(SIMILARITY_MMSEQS_PROTEIN.out.csv, expressions)

}