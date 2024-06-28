process EXTRACT_PROTEIN_SEQUENCES{

    input: 
    tuple val(meta), path(gtf), path(fasta)
    path(translated_fasta)

    output:
    tuple val(meta), path("protein_coding_gene.fa"), emit: fasta


    script: 
    """
    get_protein_coding_only.py ${gtf} protein_coding_longest_trascripts_list.txt

    # in the fasta file, modify the header and keep only transcript and gene from ...|transcript|gene|... and keep the sequence
    awk '{if(\$0 ~ /^>/){split(\$0, a, "|"); print ">"a[2]"|"a[3]}  else {print \$0}} ' ${translated_fasta} > protein_coding_modified.fa
    awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' < protein_coding_modified.fa > protein_coding_singleline.fa
    # remove th version number from the header
    sed -i 's/\\.[0-9]*//g' protein_coding_singleline.fa

    # now keep only those that are in the protein_coding_longest_trascripts_list.txt
    grep -A 1 -f protein_coding_longest_trascripts_list.txt protein_coding_singleline.fa | grep -v "^--" > protein_coding.fa

    # now remove the transcript info and keep gene only 
    awk '{if(\$0 ~ /^>/){split(\$0, a, "|"); print ">"a[2]}  else {print \$0}} ' protein_coding.fa > protein_coding_gene.fa
    """
}