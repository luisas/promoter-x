nextflow.enable.dsl = 2

include { PROMOTER_SIMILARITY_ANALYSIS } from './workflows/PROMOTER_SIMILARITY_ANALYSIS.nf'

workflow {

    // ---------------------------------------------------------------
    //              Read in the parameters
    // ---------------------------------------------------------------


    // Input parameters for PROMOTER EXTRACTION 
    gtf = Channel.fromPath("${params.gtf}").map{
        it -> [[id:"gencode"],it ]
    }
    fasta = Channel.fromPath("${params.fasta}")
    promoter_length =  params.promoter_length.toString().tokenize(',')

    gtf.combine(fasta).combine(promoter_length).multiMap{
        meta, gtf, fasta, promoter_length -> 
            files: [meta, gtf, fasta]
            promoter_length: promoter_length.toInteger()
    }.set{promoter_similarity_input}

    log.info """\
            PROMOTER SIMILARITY  ~  version 0.1"
            ======================================="
            Input sequences (FASTA)                        : ${params.fasta}
            Input gtf (annotation)                         : ${params.gtf}
            """
            .stripIndent()
    
    PROMOTER_SIMILARITY_ANALYSIS(promoter_similarity_input.files, promoter_similarity_input.promoter_length)

}