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
        "data/plier_pathways.tsv",
        "data/no_scrna_rpkm.tsv",
        "data/no_scrna_rpkm.h5",
        "data/V.tsv",
        "data/d.tsv",
        "output/Z.tsv",
        "output/U.tsv",
        "output/plier.rds",

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

# Get mouse cell type marker genes from the Marker Genes Database
rule download_marker_genes:
    output:
        "data/Mouse_cell_markers.txt"
    shell:
        "curl http://biocc.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt > data/Mouse_cell_markers.txt"

rule metadata_to_tsv:
    input:
        "data/metadata_df.rda",
        "src/1a_metadata_to_tsv.R"
    output:
        "data/recount_metadata.tsv"
    shell:
        "Rscript src/1a_metadata_to_tsv.R"

rule remove_scrna:
    input:
        "data/sra_counts.tsv",
        "data/recount_metadata.tsv",
        "src/1b_remove_scrnaseq.py"
    output:
        "data/no_scrna_counts.tsv"
    shell:
        "python src/1b_remove_scrnaseq.py "
        "data/sra_counts.tsv "
        "data/recount_metadata.tsv "
        "data/no_scrna_counts.tsv "

rule get_gene_lengths:
    input:
        "src/1_get_gene_lengths.R"
    output:
        "data/gene_lengths.tsv",
    shell:
        "Rscript src/1_get_gene_lengths.R "

rule get_pathway_matrix:
    input:
        "data/Ensembl2Reactome_All_Levels.txt",
        "data/ReactomePathwaysRelation.txt",
        "data/ReactomePathways.txt",
        "data/Mouse_cell_markers.txt",
        "src/2_create_pathway_graph.py"
    output:
        "data/plier_pathways.tsv"
    shell:
        "python src/2_create_pathway_graph.py"

rule rpkm_transform:
    input:
        "src/3_preprocess_expression.py",
        "data/no_scrna_counts.tsv",
        "data/gene_lengths.tsv",
        "data/plier_pathways.tsv"
    output:
        "data/no_scrna_rpkm.tsv"
    shell:
        "python src/3_preprocess_expression.py data/no_scrna_counts.tsv "
        "data/gene_lengths.tsv "
        "data/plier_pathways.tsv "
        "data/no_scrna_rpkm.tsv "

rule convert_to_hdf5:
    input:
        "src/4_convert_to_hdf5.R",
        "data/no_scrna_rpkm.tsv"
    output:
        "data/no_scrna_rpkm.h5"
    shell:
        "Rscript src/4_convert_to_hdf5.R"

rule calculate_pcs:
    input:
        "data/no_scrna_rpkm.tsv"
    output:
        "data/V.tsv",
        "data/d.tsv"
    shell:
        "python src/5_calculate_pcs.py data/no_scrna_rpkm.tsv data/ "

rule run_plier:
    input:
        "src/6_run_plier.R",
        "data/plier_pathways.tsv",
        "data/no_scrna_rpkm.tsv",
        "data/no_scrna_rpkm.h5",
        "data/V.tsv",
        "data/d.tsv"
    output:
        "output/plier.rds"
    shell:
        "Rscript src/6_run_plier.R"

rule save_plier_stats:
    input:
        "output/plier.rds"
    output:
        "output/Z.tsv",
        "output/U.tsv",
    shell:
        "Rscript src/7_plier_result_stats.R"