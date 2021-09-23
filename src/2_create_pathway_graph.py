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
    """
    Read the pathway file from Reactome into memory

    Arguments
    ---------
    pathway_file - The path to a tsv containing the pathway information

    Returns
    -------
    pathways - A dict that maps pathway IDs to their corresponding pathway objects
    """
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


def add_genes_to_pathways(pathways: Dict[str, Pathway], ensembl_file: str) -> Dict[str, Pathway]:
    """
    Update pathway objects by adding their corresponding genes from Reactome

    Arguments
    ---------
    pathways - A dict that maps pathway IDs to their corresponding pathway objects
    ensembl_file - The path to a tsv mapping ensembl transcript ids to pathways

    Returns
    -------
    pathways - The pathway object passed in updated to include gene information
    """
    with open(ensembl_file) as in_file:
        for line in in_file:
            line = line.strip().split('\t')

            if len(line) != 6:
                continue

            pathway_id = line[1]
            organism = line[5]
            gene_id = line[0]

            if organism != 'Mus musculus':
                continue

            if pathway_id in pathways:
                pathways[pathway_id].genes.append(gene_id)

    return pathways


def parse_tree_file(pathways: Dict[str, Pathway], tree_file: str) -> List[Tuple[str, str]]:
    """
    Read in the hierarchical relationships between pathways

    Arguments
    ---------
    pathways - A dict that maps pathway IDs to their corresponding pathway objects
    tree_file - The path to a tsv file containing parent-child pairs

    Returns
    -------
    parent_child_pairs - A list of parent - child pairs
    """
    parent_child_pairs = []
    with open(tree_file) as in_file:
        for line in in_file:
            line = line.strip().split('\t')

            # Skip empty lines or malformed lines
            if len(line) != 2:
                continue

            # Don't use non-mouse pathways
            if line[0] not in pathways or line[1] not in pathways:
                continue

            parent_child_pair = (line[0], line[1])
            parent_child_pairs.append(parent_child_pair)

    return parent_child_pairs


def build_tree(pathways: Dict[str, Pathway], tree_file: str) -> Tuple[List[TreeNode],
                                                                      Dict[str, TreeNode]]:
    """
    Assemble a tree object based on pairs of parent - child relationships

    Arguments
    ---------
    pathways - A dict that maps pathway IDs to their corresponding pathway objects
    tree_file - The path to a tsv file containing parent-child pairs

    Returns
    -------
    root_nodes - The list of pathways that have no parents
    id_to_node - A dict that points to nodes in the tree based on their pathway IDs

    """
    # Read the tree edges from the file
    child_parent_pairs = parse_tree_file(pathways, tree_file)

    root_nodes = []
    id_to_node = {}
    pair_added = set()

    for root_id in TOP_LEVEL_PATHWAYS:
        root_node = TreeNode(pathways[root_id], node_height=0)
        id_to_node[root_id] = root_node
        root_nodes.append(root_node)

    tree_modified = True

    # This is a while loop because we don't have a guarantee that the data is topo sorted
    while tree_modified:
        tree_modified = False

        for pair in child_parent_pairs:
            # Skip links we've already added
            if pair in pair_added:
                continue

            parent_id = pair[0]

            if parent_id in id_to_node:
                parent_node = id_to_node[parent_id]
                child_id = pair[1]

                # If the node already exists, add it to the list of children
                if child_id in id_to_node:
                    child_node = id_to_node[child_id]
                    if child_node.node_height > parent_node.node_height:
                        parent_node.children.append(child_node)
                # If the child node doesn't exist, make it
                else:
                    child_node = TreeNode(child_id, parent_node.node_height+1)
                    id_to_node[child_id] = child_node
                    parent_node.children.append(child_node)

                # Once the link is processed, keep track to avoid processing it again
                pair_added.add((pair))
                tree_modified = True

    return root_nodes, id_to_node


def get_leaf_nodes(nodes: List[TreeNode]) -> List[TreeNode]:
    """
    Recursively traverse pathway trees to get all the pathways without children

    Arguments
    ---------
    nodes - The nodes whose leaves you're looking for

    Returns
    -------
    leaves - The leaf nodes of the subtrees below the given nodes
    """
    leaves = []

    for node in nodes:
        # Base case
        if len(node.children) == 0:
            leaves.append(node)
        # Recursive case
        else:
            # For each node, get its leaves, and create a list appending them to each otehr
            current_leaves = get_leaf_nodes(node.children)

            leaves.extend(current_leaves)

    return leaves


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
    pathways = add_genes_to_pathways(pathways, args.ensembl_to_pathway_file)

    # Build the trees relating pathways to their children
    trees, id_to_node = build_tree(pathways, args.pathway_relation_file)

    # Select the leaf nodes
    leaves = get_leaf_nodes(trees)

    # Make sure leaves are unique
    ids_seen = set()
    unique_leaves = []
    for leaf in leaves:
        if leaf.pathway in ids_seen:
            continue
        else:
            unique_leaves.append(leaf)
            ids_seen.add(leaf.pathway)

    # TODO after finding cell type markers Create a matrix matching the PLIER format

    # Write the result