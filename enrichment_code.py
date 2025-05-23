# -*- coding: utf-8 -*-
"""Untitled1.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1_FYHDpWnNHo5eCa7mVz9Xkv0LCot10Zw
"""

!pip install gseapy pandas biomart matplotlib openpyxl

pip install mygene

import os
print(os.listdir('.'))  # Lists all files in the current directory

import os
print(os.getcwd())  # Prints the current working directory

import mygene
import gseapy as gp
import pandas as pd
import logging
from pathlib import Path

# Set up logging for error handling and debugging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Define paths to input Excel files (adjust if needed)
TURQUOISE_EXCEL_FILE = "turquoise (2).xlsx"
GREEN_EXCEL_FILE = "green_module_final (3).xlsx"
OUTPUT_DIR = "downstream_analysis_output"
Path(OUTPUT_DIR).mkdir(exist_ok=True)  # Create output directory if it doesn't exist

def load_gene_ids(turquoise_file, green_file):
    """Load Ensembl IDs from two Excel files for turquoise and green modules."""
    try:
        # Read turquoise Excel file
        turquoise_df = pd.read_excel(turquoise_file, usecols=["ENSEMBELID"])
        turquoise_ids = turquoise_df["ENSEMBELID"].dropna().tolist()
        logger.info(f"Loaded {len(turquoise_ids)} gene IDs from {turquoise_file}")

        # Read green Excel file
        green_df = pd.read_excel(green_file, usecols=["ENSEMBELID"])
        green_ids = green_df["ENSEMBELID"].dropna().tolist()
        logger.info(f"Loaded {len(green_ids)} gene IDs from {green_file}")

        return turquoise_ids, green_ids
    except FileNotFoundError as e:
        logger.error(f"File not found: {str(e)}. Please check the file paths.")
        return [], []
    except KeyError as e:
        logger.error(f"Column 'ENSEMBELID' not found: {str(e)}. Check the column name in Excel files.")
        return [], []
    except Exception as e:
        logger.error(f"Error loading Excel files: {str(e)}")
        return [], []

def annotate_genes(gene_ids, output_file, unmapped_output_file):
    """Annotate Ensembl IDs using mygene and save to CSV."""
    try:
        mg = mygene.MyGeneInfo()
        # Query mygene for annotations with returnall=True to capture unmapped
        gene_info = mg.querymany(
            gene_ids,
            scopes="ensembl.gene",
            fields=["symbol", "name", "entrezgene", "go.BP.term"],
            species="human",
            returnall=True
        )
        # Process results
        annotated_genes = []
        for g in gene_info["out"]:
            if "notfound" not in g:
                annotated_genes.append({
                    "ensembl_id": g["query"],
                    "symbol": g.get("symbol", "N/A"),
                    "name": g.get("name", "N/A"),
                    "entrezgene": g.get("entrezgene", "N/A"),
                    "go_bp": ";".join([term["term"] for term in g.get("go.BP", [])]) if g.get("go.BP") else "N/A"
                })
            else:
                annotated_genes.append({
                    "ensembl_id": g["query"],
                    "symbol": "N/A",
                    "name": "N/A",
                    "entrezgene": "N/A",
                    "go_bp": "N/A"
                })
        # Save annotated genes
        df = pd.DataFrame(annotated_genes)
        df.to_csv(output_file, index=False)
        logger.info(f"Saved annotated genes to {output_file} with {len(annotated_genes)} entries")

        # Save unmapped genes
        unmapped = gene_info["missing"]
        if unmapped:
            pd.DataFrame({"unmapped_ensembl_id": unmapped}).to_csv(unmapped_output_file, index=False)
            logger.info(f"Saved {len(unmapped)} unmapped genes to {unmapped_output_file}")
        else:
            logger.info("No unmapped genes found.")

        return df
    except Exception as e:
        logger.error(f"Error annotating genes: {str(e)}")
        return pd.DataFrame()

def perform_enrichment_analysis(symbols, module_name, output_file):
    """Perform enrichment analysis using gseapy and save results."""
    try:
        # Filter out 'N/A' symbols
        valid_symbols = [s for s in symbols if s != "N/A"]
        if not valid_symbols:
            logger.warning(f"No valid gene symbols for {module_name}. Skipping enrichment.")
            return pd.DataFrame()

        # Run_Enrichr for multiple gene sets
        enr = gp.enrichr(
            gene_list=valid_symbols,
            gene_sets=["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022", "WikiPathways_2019_Human"],
            organism="human",
            outdir=f"{OUTPUT_DIR}/{module_name}_enrichr",
            cutoff=0.5  # Relaxed p-value threshold
        )
        # Filter significant results
        results = enr.results
        significant_results = results[results["Adjusted P-value"] < 0.1]
        significant_results.to_csv(output_file, index=False)
        logger.info(f"Saved enrichment results to {output_file} with {len(significant_results)} significant terms")
        return significant_results
    except Exception as e:
        logger.error(f"Error in enrichment analysis for {module_name}: {str(e)}")
        return pd.DataFrame()

def main():
    """Main function to run gene annotation and enrichment analysis."""
    # Load gene IDs from Excel files
    turquoise_ids, green_ids = load_gene_ids(TURQUOISE_EXCEL_FILE, GREEN_EXCEL_FILE)

    if not turquoise_ids or not green_ids:
        logger.error("Failed to load gene IDs. Exiting.")
        return

    # Step 1: Annotate genes
    turquoise_annotated = annotate_genes(
        turquoise_ids,
        f"{OUTPUT_DIR}/turquoise_annotated_genes.csv",
        f"{OUTPUT_DIR}/turquoise_unmapped_genes.csv"
    )
    green_annotated = annotate_genes(
        green_ids,
        f"{OUTPUT_DIR}/green_annotated_genes.csv",
        f"{OUTPUT_DIR}/green_unmapped_genes.csv"
    )

    # Step 2: Perform enrichment analysis
    turquoise_symbols = turquoise_annotated["symbol"].dropna().tolist()
    green_symbols = green_annotated["symbol"].dropna().tolist()

    turquoise_enrichment = perform_enrichment_analysis(
        turquoise_symbols,
        "turquoise",
        f"{OUTPUT_DIR}/turquoise_enrichment.csv"
    )
    green_enrichment = perform_enrichment_analysis(
        green_symbols,
        "green",
        f"{OUTPUT_DIR}/green_enrichment.csv"
    )

    # Check success
    if not turquoise_annotated.empty and not green_annotated.empty and \
       not turquoise_enrichment.empty and not green_enrichment.empty:
        logger.info("Downstream analysis completed successfully!")
    else:
        logger.warning("Some analyses failed. Check logs for details.")

if __name__ == "__main__":
    main()