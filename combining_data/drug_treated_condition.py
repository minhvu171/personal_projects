import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# FUNCTION TO LOAD AND CLEAN A DATASET
def load_and_clean(filepath, sheet_name=0, skip_rows=0, gene_col=2, logfc_col=4, padj_col=5):
    """
    Reads an Excel sheet, extracts specific columns, removes rows with NA gene names.

    Parameters:
        filepath (str): Path to the Excel file
        sheet_name (int or str): Sheet index or name to read
        skip_rows (int): Number of rows to skip at the top
        gene_col, logfc_col, padj_col (int): Column indices (0-based) to extract

    Returns:
        pd.DataFrame: Cleaned dataset with columns: Gene, Log2FC, AdjPValue
    """
    # Read the sheet
    df = pd.read_excel(filepath, sheet_name=sheet_name, skiprows=skip_rows)

    # Select required columns by index
    df = df.iloc[:, [gene_col, logfc_col, padj_col]]

    # Rename columns for consistency
    df.columns = ["Gene", "Log2FC", "AdjPValue"]

    # Remove rows where Gene is NA
    df = df[df["Gene"].notna() & (df["Gene"] != "NA")]

    # Reset index
    df.reset_index(drop=True, inplace=True)

    return df


# LOAD FIRST TWO DATASETS (Unfortunately, these are all data from my university, so I can't share it)
# The R234 gene file is structured differently from the first two, so I had 
# to process it seperately.
DBT = "DBT.xlsx"
drug5A = "5A.xlsx"
R234 = "R234.xlsx"

DBT_df = load_and_clean(DBT, sheet_name=0, skip_rows=1)  
drug5A_df = load_and_clean(drug5A, sheet_name=0, skip_rows=0)  

print(f"DBT File : {len(DBT_df)} genes")
print(f"5A File : {len(drug5A_df)} genes")

# LOAD AND MERGE SHEETS IN R234 DATASET 
R234_df_sheet3 = load_and_clean(R234, sheet_name=2, skip_rows=1, gene_col=2, logfc_col=5, padj_col=6)
R234_df_sheet4 = load_and_clean(R234, sheet_name=3, skip_rows=1, gene_col=2, logfc_col=5, padj_col=6)

# Combine both sheets
R234_df = pd.concat([R234_df_sheet3, R234_df_sheet4], ignore_index=True)
print(f"R234: {len(R234_df)} genes") # Checking whether reading data is succesful


# SAVE CLEANED DATASETS
#DBT_df.to_excel("DBT_cleaned.xlsx", index=False)
#drug5A_df.to_excel("5A_cleaned.xlsx", index=False)
#R234_df.to_excel("R234_cleaned.xlsx", index=False)

# FIND COMMON GENES AMONG THE THREE LISTS OF GENES
DBT_genes = set(DBT_df["Gene"])
drug5A_genes = set(drug5A_df["Gene"])
R234_genes = set(R234_df["Gene"])

common_genes = DBT_genes & drug5A_genes & R234_genes
print(f"Common genes: {len(common_genes)}")


# CREATE VENN DIAGRAM
plt.figure(figsize=(6,6))
venn3([DBT_genes, drug5A_genes, R234_genes],
      set_labels=("DBT", "5A", "R234"))
plt.title("Gene Overlap between the three drug-treated conditions")
plt.savefig("drug_gene_overlap_venn.png")
plt.close()

# SAVE COMMON GENES WITH ALL VALUES
# Merge on "Gene" to keep all related values
common_df = pd.DataFrame({"Gene": list(common_genes)})
common_df = common_df.merge(DBT_df, on="Gene", how="left")
common_df = common_df.merge(drug5A_df, on="Gene", how="left", suffixes=("", "_5A"))
common_df = common_df.merge(R234_df, on="Gene", how="left", suffixes=("", "_R234"))
#common_df.to_excel("Cell_lines_common_genes.xlsx", index=False)

