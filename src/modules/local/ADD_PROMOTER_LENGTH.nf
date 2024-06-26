process ADD_PROMOTER_LENGTH{

    input: 
    tuple val(meta), path(csv)
    val(promoter_length)

    output:
    tuple val(meta), path("sim_length.csv"), emit: csv
    
    script:
    """
    # add promoter length to the csv file as a new column
    # first line is the header, append the string "promoter_length" to it
    # then append the promoter length to each line
    head -n 1 $csv | awk -v OFS=',' '{print \$0, "promoter_length"}' > sim_length.csv
    tail -n +2 $csv | awk -v promoter_length=${promoter_length} 'BEGIN{FS=OFS=","} {print \$0, promoter_length}' >> sim_length.csv
    """

}