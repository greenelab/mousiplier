rule all:
    input:
        "data/sra_counts.tsv",
        "data/metadata_df.rda",
        "data/recount_metadata.tsv",
        "data/no_scrna_counts.tsv",
        "data/gene_lengths.tsv",
        "data/Ensembl2Reactome_All_Levels.txt",
        "data/ReactomePathwaysRelation.txt",
        "data/ReactomePathways.txt",

rule download_data:
    output:
        "data/sra_counts.tsv",
        "data/recount_metadata.rda"
    shell:
        "Rscript src/0_download_recount3.R "

rule download_reactome_pathways:
    output:
        "data/Ensembl2Reactome_All_Levels.txt",
        "data/ReactomePathwaysRelation.txt",
        "data/ReactomePathways.txt"
    shell:
        "curl https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt > data/Ensembl2Reactome_All_Levels.txt ; "
        "curl https://reactome.org/download/current/ReactomePathways.txt > data/ReactomePathways.txt ; "
        "curl https://reactome.org/download/current/ReactomePathwaysRelation.txt > data/ReactomePathwaysRelation.txt "

rule metadata_to_tsv:
    input:
        "data/metadata_df.rda"
    output:
        "data/recount_metadata.tsv"
    shell:
        "Rscript src/metadata_to_tsv.R"

rule remove_scrna:
    input:
        "data/sra_counts.tsv",
        "data/recount_metadata.tsv"
    output:
        "data/no_scrna_counts.tsv"
    shell:
        "python src/remove_scrnaseq.py "
        "data/sra_counts.tsv "
        "data/recount_metadata.tsv "
        "data/no_scrna_counts.tsv "

rule get_gene_lengths:
    output:
        "data/gene_lengths.tsv"
    shell:
        "Rscript src/1_get_gene_lengths.R "