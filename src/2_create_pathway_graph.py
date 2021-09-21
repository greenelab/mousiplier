"""
This script parses a number of files from Reactome to match pathways to their
corresponding titles and genes
"""

import argparse
from dataclasses import dataclass, field
from typing import List, Dict, Tuple

# Top level pathways from https://reactome.org/PathwayBrowser as of 9/21
TOP_LEVEL_PATHWAYS = ['R-MMU-9612973', 'R-MMU-1640170', 'R-MMU-1500931', 'R-MMU-8953897',
                      'R-MMU-4839726', 'R-MMU-1266738', 'R-MMU-8963743', 'R-MMU-73894',
                      'R-MMU-69306', 'R-MMU-1474244', 'R-MMU-74160', 'R-MMU-109582',
                      'R-MMU-168256', 'R-MMU-1430728', 'R-MMU-392499', 'R-MMU-8953854',
                      'R-MMU-397014', 'R-MMU-112316', 'R-MMU-1852241', 'R-MMU-5357801',
                      'R-MMU-9609507', 'R-MMU-1474165', 'R-MMU-9709957', 'R-MMU-162582',
                      'R-MMU-382551', 'R-MMU-5653656'
                      ]

@dataclass
class Pathway():
    id: str
    name: str
    genes: List[str] = field(default_factory=list)


@dataclass
class TreeNode():
    pathway: Pathway
    # Need to keep track of node height to make sure trees don't grow cycles
    node_height: int
    children: List[Pathway] = field(default_factory=list)


def parse_pathway_file(pathway_file:str) -> Dict[str, Pathway]:
    pathways = {}
    with open(pathway_file) as in_file:
        for line in in_file:
            line = line.strip().split('\t')
            id = line[0]
            name = line[1]
            organism = line[2]

            # We only want mouse pathways
            if organism != 'Mus musculus':
                continue

            current_pathway = Pathway(id, name)

            pathways[id] = current_pathway

    return pathways


def add_genes_to_pathways(pathway: Dict[str, Pathway], ensembl_file: str) -> Dict[str, Pathway]:
    pass


def parse_tree_file(tree_file: str, pathways: Dict[str, Pathway]) -> List[Tuple[str, str]]:
    parent_child_pairs = []
    with open(tree_file) as in_file:
        for line in in_file:
            line = line.strip().split('\t')

            # Skip empty lines or malformed lines
            if len(line != 2):
                continue

            # Don't use non-mouse pathways
            if line[0] not in pathways or line[1] not in pathways:
                continue

            parent_child_pair = (line[0], line[1])
            parent_child_pairs.append(parent_child_pair)

    return parent_child_pairs


def build_tree(pathways: Dict[str, Pathway], tree_file: str) -> Tuple(List[TreeNode],
                                                                      Dict[str: TreeNode]):
    """Returns the root nodes and a dict mapping ids to pointers within the tree"""
    child_parent_pairs = parse_tree_file(tree_file)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pathway_file',
                        help='The ReactomePathways.txt file listing each pathway '
                             'and its corrsesponding id and species',
                        default='data/ReactomePathways.txt')
    parser.add_argument('--pathway_relation_file',
                        help='The ReactomePathwaysRelation.txt that lists each '
                             'pathway - child pair',
                        default='data/ReactomePathwaysRelation.txt')
    parser.add_argument('--pathway_gene_file',
                        help='The file mapping each pathway to its corresponding genes',
                        default='data/pathway_to_gene.tsv')
    parser.add_argument('--ensembl_to_pathway_file',
                        help='A file mapping ensembl ids to reactome pathways',
                        default='data/Ensembl2Reactome_All_Levels.txt')
    parser.add_argument('--out_file',
                        help='The file to store the selected pathways for use in PLIER',
                        default='data/plier_pathways.tsv')
    args = parser.parse_args()

    # Read all pathways, store the mouse pathways
    pathways = parse_pathway_file(args.pathway_file)

    # Find which genes are in each pathway
    pathways = add_genes_to_pathways(args.ensembl_to_pathway_file)

    # Build the trees relating pathways to their children
    trees, id_to_node = build_tree(pathways)

    # Find which pathways to keep

    # Create a matrix matching the PLIER format

    # Write the result