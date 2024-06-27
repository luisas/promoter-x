process ADD_EXPRESSION{

    input: 
    tuple val(meta), path(csv)
    tuple val(meta), path(expression_files)

    output:
    tuple val(meta), path("${meta.id}_exp.csv"), emit: csv

    script: 
    """
    02_add_expression.py $csv $expression_files "${meta.id}_exp.csv"
    """
}