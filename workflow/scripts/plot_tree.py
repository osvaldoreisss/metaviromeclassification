import argparse
from Bio import Phylo
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tree_file')
    ap.add_argument('--out_file')

    args = ap.parse_args()

    tree = Phylo.read(args.tree_file, "newick")

    fig = plt.figure(figsize=(10, 60), dpi=100)
    axes = fig.add_subplot(1, 1, 1)
    
    Phylo.draw(tree, axes=axes, do_show=False)

    plt.savefig(args.out_file)

if __name__ == "__main__":
    main()