"""
Runs the parameters for a given simulation in an input
pickle file and writes the output tree sequence to
an output file.
"""
import pickle
import sys

import msprime


if __name__ == "__main__":
    assert msprime.__version__ == "0.7.4"

    if len(sys.argv) != 4:
        raise ValueError(
            "Usage: python3 run_old_msprime.py <num_replicates> "
            "<params file> <output prefix>"
        )

    num_replicates = int(sys.argv[1])
    with open(sys.argv[2], "rb") as f:
        kwargs = pickle.load(f)
    output_prefix = sys.argv[3]

    kwargs["num_replicates"] = num_replicates
    for j, ts in enumerate(msprime.simulate(**kwargs)):
        ts.dump(f"{output_prefix}_{j}.trees")
