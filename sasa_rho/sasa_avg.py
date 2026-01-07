import numpy as np
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description="Average SASA")
parser.add_argument("-i", type=str, help="input path")
parser.add_argument("-o", type=str, help="ouput filename")
args = parser.parse_args()

files           = [args.i]
output_filename = args.o

data = defaultdict(list)
for fname in files:
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # skip header / comment lines
            if line.startswith("Computed") or line.startswith("#"):
                continue

            res, value = line.split()
            data[res].append(float(value))

# write output
with open(output_filename, "w") as out:
    out.write("Residue  Mean  StdError\n")
    for res in sorted(data):
        values = np.asarray(data[res], dtype=float)
        n = len(values)

        if n == 0:
            mean = np.nan
            stderr = np.nan

        elif n == 1:
            mean = values[0]
            stderr = 0.0   # no uncertainty with one sample

        else:
            mean = values.mean()
            stderr = values.std(ddof=1) / np.sqrt(n)

        out.write(f"{res:4s} {mean:8.3f} {stderr:8.3f}\n")

