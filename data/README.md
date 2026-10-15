# Proteomic Datasets

The raw SWATH-MS proteomic matrices and statistical analysis tables used in this project are available from the author upon reasonable request.

To reproduce the analysis from scratch, place the corresponding files into this directory:
- `swathms_mutants.csv`: Filtered SWATH-MS protein abundance matrix for *S. cerevisiae* mutant and wild-type strains.
- `final_merged_data.csv`: Merged and codon-annotated dataset across all conditions.
- `limma/`: Directory containing differential expression results (*limma* summary tables for the 5 biological contrasts).
- `kegg/`: Directory containing KEGG pathway enrichment output tables.

*Note: All thesis figures are already pre-computed and viewable directly in `Minimalist_Plot_Generator.ipynb` without needing the raw data files.*
