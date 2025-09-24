from utils.unzip import unzip_data
from utils.get_genes import create_table_with_genes

from base_analysis import density_plot
from diff_analysis import differential_analysis

COMPRESSED_DATA_PATH = "src/data/SRP120552.zip"
METADATA_PATH = "src/data/metadata_SRP120552.tsv"

def main() -> None:
    # Unzip data and store in data folder
    data_path = unzip_data(COMPRESSED_DATA_PATH)
    
    # --- PART 1 ---
    # b. Get gene names
    print("Getting gene names and mapping...")
    data_path = create_table_with_genes(data_path)

    # c. Meta data analysis
    differential_analysis(data_path, METADATA_PATH)
    # print("Running metadata analysis on table...")
    # density_plot(data_path)

    # data had a lower variance of
    # about 0.85 with the total range of log scaled data spanning from about 0.45 to 7.58. A vast majority of the
    # data is centralized between 0 and 1 (about 80%) with some outliers skewing the data towards the right.
    

if __name__ == "__main__":
    main()