import numpy as np
from collections import defaultdict

files = ["1XQ8_SASA","2NBI_SASA","6XRY_SASA"]

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
with open("averaged_sasa", "w") as out:
    out.write("Residue  Mean  StdError\n")
    for res in sorted(data):
        values = np.array(data[res])
        mean = values.mean()
        stderr = values.std(ddof=1) / np.sqrt(len(values))
        out.write(f"{res:4s} {mean:8.3f} {stderr:8.3f}\n")

