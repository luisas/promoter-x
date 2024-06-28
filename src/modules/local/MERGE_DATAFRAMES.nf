process  MERGE_DATAFRAMES{

    input:
    tuple val(meta), path(csv1)
    tuple val(meta2), path(csv2)

    output:
    tuple val(meta), path("${meta.id}_merged.csv"), emit: csv

    script:
    """
    merge_pair_files.py $csv1 $csv2 "${meta.id}_merged.csv"
    """


}