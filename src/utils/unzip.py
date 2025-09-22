import gzip
from pathlib import Path
import shutil

def unzipData(source):

    path = Path(source)
    if not path.with_suffix('').exists():
        dest = path.with_suffix('')
        with gzip.open(path, 'rb') as f_in, open(dest, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

