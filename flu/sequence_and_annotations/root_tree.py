'''
root_tree.py

takes a newick tree file as argument, loads it as a biopython tree
if second argument is provided, looks for a leaf that starts with this argument
otherwise tries to root with A/Hong_Kong/68
'''
from Bio import Phylo
import sys

if len(sys.argv)>1:
    tree_name = sys.argv[1]
else:
    tree_name = '../data/H3N2_HA1_all_years_filtered.nwk'

if len(sys.argv)>2:
    outgroup_mask = sys.argv[2]
else:
    outgroup_mask = 'A/Hong_Kong/68'


HA_tree = Phylo.read(tree_name, 'newick')
for leaf in HA_tree.get_terminals():
    leaf.name = '/'.join(leaf.name.split('/')[:-1])

HA_outgroup = [leaf for leaf in HA_tree.get_terminals() if leaf.name.startswith(outgroup_mask)][0]
HA_tree.root_with_outgroup(HA_outgroup)
HA_tree.ladderize()

Phylo.write(HA_tree, tree_name, 'newick')

