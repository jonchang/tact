#!/usr/bin/env python3
"""Script to download and set up the Percomorphaceae example dataset."""

import csv
import os
import subprocess
import sys
import tempfile
import urllib.request

EXAMPLES_DIR = os.path.dirname(os.path.abspath(__file__))


def download_file(url, dest):
    """Download a file from URL to destination."""
    print(f"Downloading {os.path.basename(dest)}...")
    urllib.request.urlretrieve(url, dest)


def main():
    """Sets up the dataset."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Download phylogeny
        phylogeny_url = "https://fishtreeoflife.org/downloads/taxonomy/subdivision/Percomorphaceae.tre"
        phylogeny_path = os.path.join(EXAMPLES_DIR, "Percomorphaceae.tre")
        download_file(phylogeny_url, phylogeny_path)

        # Download and decompress taxonomy
        taxonomy_url = "https://fishtreeoflife.org/downloads/PFC_taxonomy.csv.xz"
        taxonomy_xz = os.path.join(temp_dir, "PFC_taxonomy.csv.xz")
        download_file(taxonomy_url, taxonomy_xz)

        print("Decompressing taxonomy file...")
        subprocess.run(["unxz", taxonomy_xz], check=True)
        taxonomy_csv = os.path.join(temp_dir, "PFC_taxonomy.csv")

        # Filter to Percomorphaceae entries
        print("Filtering taxonomy to Percomorphaceae entries...")
        filtered_rows = []
        with open(taxonomy_csv, encoding="utf-8") as f:
            reader = csv.reader(f)
            header = next(reader)
            filtered_rows.append(header)

            # Find indices for columns to remove
            try:
                superclass_idx = header.index("superclass")
                division_idx = header.index("division")
            except ValueError:
                print("Warning: Could not find superclass or division columns", file=sys.stderr)
                superclass_idx = None
                division_idx = None

            # Filter rows containing Percomorphaceae
            for row in reader:
                if any("Percomorphaceae" in cell for cell in row):
                    filtered_rows.append(row)

        # Remove columns from superclass to division (inclusive)
        print("Removing unneeded ranks and filling empty cells...")
        if superclass_idx is not None and division_idx is not None:
            cols_to_remove = set(range(superclass_idx, division_idx + 1))
        else:
            cols_to_remove = set()

        # Process rows: remove columns and fill empty cells
        processed_rows = []
        for row_idx, row in enumerate(filtered_rows):
            if row_idx == 0:
                # Header row
                new_header = [col for i, col in enumerate(row) if i not in cols_to_remove]
                processed_rows.append(new_header)
                col_names = new_header
            else:
                # Data row
                new_row = [cell for i, cell in enumerate(row) if i not in cols_to_remove]

                # Fill empty cells with value from right + column name
                for i in range(len(new_row) - 1):  # Don't process last column (species name)
                    if not new_row[i].strip():  # Empty cell
                        # Find next non-empty cell to the right
                        for j in range(i + 1, len(new_row)):
                            if new_row[j].strip():  # Found non-empty cell
                                new_row[i] = f"{new_row[j]}_{col_names[i]}"
                                break

                processed_rows.append(new_row)

        # Write output
        output_path = os.path.join(EXAMPLES_DIR, "Percomorphaceae.csv")
        with open(output_path, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(processed_rows)

        print("Setup complete! Files saved to:")
        print(f"  - {phylogeny_path}")
        print(f"  - {output_path}")


if __name__ == "__main__":
    main()
