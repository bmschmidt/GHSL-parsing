from quadfeather.tiler import Quadtree, Rectangle
from pathlib import Path
import duckdb
import typer
from pyarrow import parquet as pq
import numpy as np
from numpy import random
from pyarrow import feather
import pyarrow as pa

app = typer.Typer()


@app.command()
def shuffle(dest: Path):
    if dest.with_suffix(".feather").exists():
        return
    print(f"Shuffling {dest}")
    ys = [(str(y), "bool") for y in range(1975, 2027, 5)]
    r = pq.read_table(dest, schema=pa.schema([("x", "float32"), ("y", "float32"), *ys]))

    indices = np.arange(r.num_rows)
    random.shuffle(indices)
    r = r.take(indices)
    feather.write_feather(r, dest.with_suffix(".feather"))


if __name__ == "__main__":
    app()
