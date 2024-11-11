from quadfeather.tiler import Quadtree, Rectangle
from pathlib import Path
import typer
import pyarrow as pa
from pyarrow import ipc, compute as pc
import random

app = typer.Typer()


#
def batch_sizes(p):
    if p.with_suffix(".batch_sizes.txt").exists():
        return [*map(int, p.with_suffix(".batch_sizes.txt").open().readlines())]
    r = ipc.open_file(p)
    with p.with_suffix(".batch_sizes.txt").open("w") as fout:
        for i in range(r.num_record_batches):
            batch = r.get_batch(i)
            fout.write(f"{batch.num_rows}\n")
    return [*map(int, p.with_suffix(".batch_sizes.txt").open().readlines())]


@app.command()
def build(dest: Path):

    print(f"Building quadtree at {dest}")
    tree = Quadtree(
        basedir=dest.expanduser(),
        extent=Rectangle((-20040000, 20040000), (-14364430, 7536672)),
        mode="write",
    )
    total_slices = 100
    slice_num = 0

    while slice_num < total_slices:
        local_batch = []
        for p in Path("points").glob("*.feather"):
            print(f"Reading {p}")
            counts = batch_sizes(p)
            if len(counts) == 0:
                continue
            tb = ipc.open_file(p)

            total_in_batch = sum(counts)
            slice_size = int(total_in_batch / total_slices)
            start_ix = slice_size * slice_num
            end_ix = (slice_num + 1) * slice_size

            batch_start = 0
            for i, batch_length in enumerate(counts):

                if batch_start > end_ix:
                    break
                batch_end = batch_start + batch_length

                # Calculate the overlap between the desired indices and the current batch
                start_overlap = max(start_ix, batch_start)
                end_overlap = min(end_ix, batch_end)
                if start_overlap >= end_overlap:
                    batch_start = batch_end
                    continue

                batch = tb.get_batch(i)
                start_within_batch = start_overlap - batch_start
                end_within_batch = end_overlap - batch_start

                # Generate chunk boundaries within the batch
                chunk_size = 32
                chunk_starts = list(
                    range(start_within_batch, end_within_batch, chunk_size)
                )
                chunk_ends = chunk_starts[1:] + [end_within_batch]

                for a, b in zip(chunk_starts, chunk_ends):
                    subset = batch.slice(a, b - a)
                    local_batch.append(subset)

                batch_start = batch_end

        print(f"Number of local batches: {len(local_batch)}")

        random.shuffle(local_batch)
        if len(local_batch) == 0:
            raise ValueError("No data")

        tb_batches = [[]]
        current_length = 0
        for batch in local_batch:

            if current_length + batch.num_rows > 65536:
                tb_batches.append([])
                current_length = 0
            tb_batches[-1].append(batch)
            current_length += batch.num_rows

        tb = pa.concat_tables(
            [pa.Table.from_batches(b).combine_chunks() for b in tb_batches]
        )

        # I forgot that we use browser-style coordinates so need to flip here.
        tb = tb.drop(["y"]).append_column("y", pc.multiply(tb["y"], -1))
        print("Inserting data into quadtree for " + str(slice_num))
        print(
            f"Number of rows: {tb.num_rows}, number of batches: {len(tb.to_batches())}, average length of batch: {tb.num_rows / len(tb.to_batches())}"
        )
        tree.insert(tb)
        print("Data insertion complete" + str(slice_num))
        slice_num += 1
    print("Finalizing quadtree")
    tree.finalize()
    print("Quadtree build complete")


if __name__ == "__main__":
    app()
