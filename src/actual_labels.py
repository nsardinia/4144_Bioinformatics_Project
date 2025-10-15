from pathlib import Path

import pandas as pd
import numpy as np

METADATA_PATH = Path("data\metadata_SRP120552.tsv").absolute()

labels_list = ['t0', 't24', 't48', 't72']

def get_actual_labels() -> list[int]:
    meta = pd.read_csv(METADATA_PATH, sep="\t")

    metadata = pd.DataFrame(
        zip(
            meta["refinebio_accession_code"], 
            meta["refinebio_time"],
        ),
        columns=["sample", "condition"]
    )

    return metadata["condition"].map(lambda x: labels_list.index(x)).to_numpy()