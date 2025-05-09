process ADD_EXPRESSION{

    input: 
    tuple val(meta), path(csv)
    tuple val(meta), path(expression_files)

    output:
    tuple val(meta), path("${meta.id}_exp.csv"), emit: csv
    path("${meta.id}_samples_mean.csv"), emit: samples_mean

    script: 
    """
    02_add_expression_vector.py $csv $expression_files "${meta.id}_exp.csv"
    compute_sample_geometric_mean.py $expression_files "${meta.id}_samples_mean.csv"
    """
}