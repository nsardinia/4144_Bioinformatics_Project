import zipfile
from pathlib import Path
import shutil


def unzip_data(source: str) -> str:
    """Unzip the data file."""
    path = Path(source).absolute()
    
    if not path.with_suffix('').exists():
        dest = path.with_suffix('')
        
        with zipfile.ZipFile(path, 'r') as zip_ref:
            zip_ref.extract("SRP120552.tsv", path=dest.parent)
        
        return dest.parent / "SRP120552.tsv"

    return path           
