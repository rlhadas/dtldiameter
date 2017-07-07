# NewDiameter.py
# Written by Eli Zupke and Andrew Ramirez July 2017
# It is based off of Jordan Haack's work on a polynomial-time algorithm for the DTL MPR Diameter Problem

# TODO: Update header

# TODO: provide ancestry details.

# 0. ON THE PAPER
#
#       This program is based off of an algorithm written by Jordan Haack. There is a paper that describes this
#   algorithm more thoroughly than can be expressed in source code comments. The paper title is "Computing the Diameter
#   of the Space of Maximum Parsimony Reconciliations in the Duplication-Transfer-Loss Model". #TODO add more info
#
#       Note that in this documentation and in that paper, there is a notion of the "group" of a gene node that contains
#   every mapping and event node relating to that gene. While this is a useful abstraction to think about the problem,
#   it doesn't seem to be required in the code. So when the documentation of this program refers to the "group" of a
#   node, sometimes notated as Group(u), be aware that it does not directly correspond to any data structure in the
#   program.

# 1. ON TREE REPRESENTATION FORMATS:
#
#       This file deals with trees in two formats: Edge-based formats, and vertex-based formats. The edged-based
#   trees are what are output from DTLReconGraph.py (which returns them straight from newickFormatReader), and are what this
#   file uses to represent the species tree. For readability and convenience purposes, the gene tree is converted to
#   a vertex tree in the reformat_tree() method.
#
#   Both trees use dictionaries to store the tree structure, but formats differ slightly.
# Example:
#   When N is a node descended from R, with children C1 and C2,
#
#   An entry in the edge-based representation looks like this:
#       {('R','N'): ('R','N', ('N','C1'), ('N','C2')) ...}
#
#   And an entry in the vertex-based representation looks like this:
#       {'N':('C1','C2') ...}
#
#   Because the standard elsewhere in this codebase is to use edge-based trees, all vertex-based tree variable names
#   will explicitly contain the word "vertex", while edge-based trees will not. So, an edge-based tree can be named
#   "tree" while a vertex-based tree has to be named "vertex_tree".

# 2. ON THE NAMING CONVENTION OF THE TWO TREES:
#
#       This file calls the two trees the gene tree and the species tree. Other programs, like DTLReconGraph.py, have different
#   naming conventions (such as "host" for the species tree and "parasite" for the gene tree) because they were coded
#   under different assumptions as to what the two trees represent. Let it be understood that these names are
#   synonymous, and that references to "hTop" or "pTop" refer to the name of the handle of the species and
#   gene trees that DTLReconGraph outputs.


# 3. ON THE USAGE OF PATHS AND THE PATH REPRESENTATION FORMAT:
#
#       compute_enter_mapping_table() needs to use the symmetric set difference between every pair of paths on the gene
#   tree. This is pre-computed in build_path_symmetric_set_difference_table() and stored in a nested dictionary.
#   The format of each path is a tuple where the first entry is the species node that the path starts on, and the second
#   entry is the species node that the path ends on. These paths are stored as the keys of a dictionary called
#   path_edges.

# 4. ON THE DYNAMIC PROGRAMMING TABLES AND THEIR FUNCTIONS:
#
#       The Diameter algorithm involves the use of three dynamic programming tables:
#   The enter_mapping_table, the exit_mapping_table, and the exit_event_table. All three of these have the same key
#   format: [u][x][y], where u is a node on the gene tree, and x and y are either event nodes or mapping nodes. Which
#   type each table has is evident in the name.
#
#       These tables are computed for every gene node, in post-order. Functions that begin with the word 'computed'
#   build these tables by taking them in as a reference and modifying the values. These functions return nothing.
#   Other functions (such as one starting with 'build') actually return the thing that they are creating.
#
#       The tables are described in detail in the paper, but a summary of their functions follows:
#
#       exit_event_table (only implied in the paper) contains the largest number of event nodes that each pair of
#   reconciliations rooted at events x and y can differ for, where x and y are exit-event nodes in Group(u).
#       exit_mapping_table (called EXIT in the paper) contains the largest number of event nodes that each pair of
#   reconciliation subtrees rooted at mapping nodes x and y can differ for, on the condition that x and y go
#   immediately to an exit-event.
#       enter_mapping_table (called ENTER in the paper) contains the largest number of event nodes that each pair
#   of reconciliation subtrees rooted at x and y can differ for when x and y are used to enter Group(u).
#
#       To get a better sense of how the tables are structured, you can try running an example in debug mode, which
#   will actually print out these tables.

# 5. ON EXIT_X_BY_Y DICTIONARIES:
#
#       There are several helper dictionaries which pertain to exit events. Their names should be pretty self
#   explanatory, but they are described more fully here.
#
#       exit_events_by_mapping is keyed by mapping node. Each entry contains a list of every exit event
#   that is a child of that mapping node
#       exit_events_by_gene is keyed by gene node. Each entry contains a list of every exit event in
#   the group of that gene node
#       exit_mappings_by_gene is keyed by gene node. Each entry contains a list of every mapping node
#   in group of that gene node, on the condition that that mapping node has an exit event

# -1. DATA STRUCTURE QUICK REFERENCE:
#
#   Ancestor Comparability:
#       A ancestor of B:    'an'
#       A descendant of B:  'des'
#       A incomparable to B:'in'
#       A is equal to B:    'eq'
#
#   Pre clean_graph():
#
#      DTL Reconciliation graph:
#           { mapping_node: [event1, event2, ... eventn, number] ...}
#
#      Event node:
#           ('type', child_mapping_node1, child_mapping_node2, number)
#
#   Post clean_graph():
#
#       DTL Reconciliation graph (post clean_graph):
#           { mapping_node: [event1, event2, ...] ...}
#
#       Event node (post clean_graph):
#           ('type', child_mapping_node1, child_mapping_node2)
#
#
#   Mapping node:
#       ('gene_node','SPECIES_NODE')
#   or in loss or contemporary event nodes:
#       (None, None)
#
#
#   (edge) trees:
#       {('R','N'): ('R','N', ('N','C1'), ('N','C2')) ...}
#       aka:
#       {root_edge: (root_edge[0], root_edge[1], child1_edge, child2_edge) ...}
#
#   vertex_trees:
#       {'N':('C1','C2') ...}
#
#   path:
#       ('SOURCE','DESTINATION')
#
#   path_edges:
#       {path: [edge1, edge2, edge3, ...] ...}
#
#   exit_events_by_mapping:
#       {mapping_node: [exit_event1, exit_event2, ...] ...}
#
#   exit_events_by_gene:
#       {gene_node: [exit_event1, exit_event2, ...] ...}
#
#   exit_mappings_by_gene:
#       {gene_node: [exit_mapping_node1, exit_mapping_node2, ...] ...}

import DTLReconGraph
import time
import csv
import os.path
from collections import OrderedDict
from itertools import product
# Used for command line arguments:

import sys
import re
import traceback
import optparse


def reformat_tree(tree, root):
    """A recursive function that changes the format of a (species or gene) tree from edge to vertex, as described
    above. It returns the tree (in postorder), the root of the tree, and the number of nodes in the tree. The base
    case of this function is when there are no children for a given root.
    :param tree:        A tree in edge format
    :param root:        The root of that tree
    :return:            0: The new vertex based tree,
                        1: The root of that tree, and
                        2: The number of nodes in that tree. """

    # This line catches the "xTop" handle and replaces
    new_root = root[1] if isinstance(root, tuple) else tree[root][1]

    child1 = tree[root][2][1] if tree[root][2] is not None else None  # These lines handle the leaves, where
    child2 = tree[root][3][1] if tree[root][3] is not None else None  # there is None in the place of a tuple

    # This is the tree that we will be returning. We will add the subtrees of our children first, then add this node.
    new_vertex_tree = OrderedDict()  # This has to be an OrderedDict, otherwise we can't guarantee it's in postorder

    # This is the number of nodes in the subtree rooted at each of our children
    # We actually don't need to calculate this here at all, since the count will just be the length of new_vertex_tree
    child1_count = 0
    child2_count = 0

    if child1 is not None:  # If this node has children, then we need to add their children's subtrees to the dict
        child1_tree, _, child1_count = reformat_tree(tree, tree[root][2])
        new_vertex_tree.update(child1_tree)
    if child2 is not None:
        child2_tree, _, child2_count = reformat_tree(tree, tree[root][3])
        new_vertex_tree.update(child2_tree)
    new_vertex_tree.update({new_root: (child1, child2)})  # We add this node last, to ensure postorderness.

    return new_vertex_tree, new_root, (child1_count + child2_count + 1)

def find_paths_ending_past_here(current_edge, previous_edges, tree):
    """A recursive function to find all of the valid paths through a binary edge-based tree that end at or below a
    currently examined edge, and the edges included in each path.
     :param current_edge:   The edge we are examining in this function call
     :param previous_edges: A list of the edges that you would traverse to get to current_edge from the root of the
                            tree, in order
     :param tree:           The (edge-based) tree in its entirety.
     :return:               A dictionary that is keyed by path and contains lists of edges for each path."""

    # This list contains every node we visited to get here from the current_edge
    next_edges = previous_edges + [current_edge]
    path_edges = {}  # This is the edges contained in each path.
    for i in range(0, len(next_edges)):
        source_node = next_edges[i]
        # new_path becomes every path that ends in this value, including A->A
        if i >= len(next_edges) - 1:
            new_path = (tree[source_node][1], tree[source_node][1])  # An ugly hack to account for the top of the tree
        else:
            new_path = (tree[source_node][1], current_edge[1])  # This is what we should always be able to do
        path_edges[new_path] = next_edges[i + 1:len(next_edges)]  # This is the list of edges

    child1 = tree[current_edge][2]
    child2 = tree[current_edge][3]

    if child1 is not None:  # Then this Node is not a leaf Node, so we need to add this Node's children
        child1_path_edges = find_paths_ending_past_here(child1, next_edges, tree)
        child2_path_edges = find_paths_ending_past_here(child2, next_edges, tree)
        child1_path_edges.update(child2_path_edges)
        path_edges.update(child1_path_edges)
        # Otherwise, we have reached the end of the tree (the base case)
    return path_edges


def find_valid_paths(root, tree):
    """A function that uses find_paths_ending_past_here to find all of the valid paths on the given edge-based tree
    and returns the list of paths, a dict containing the edges for each path, and a list of
    non-trivial paths (where source != destination)
    :param root:        The root node of the tree
    :param tree:        The (edge-based) tree to find valid paths through
    :return:            0: A dict containing the lists of edges contained in each vaild path across the tree, and
                        1: A list of non-trivial paths (where the source is not the same as the destination) """

    path_edges = find_paths_ending_past_here(root, [], tree)

    # Strip out any trivial paths where source = destination (such as A->A). They have no path edges, so the length is
    # zero.
    non_trivial_path_list = filter(lambda x: len(path_edges[x]) > 0, path_edges.keys())
    return path_edges, non_trivial_path_list


def build_path_symmetric_set_difference_table(path_edges):
    """Computes the table containing the number of nodes in the symmetric set difference between any two paths on the
     species tree sTree. This is used in assigning a score to two paths' losses.
     :param path_edges:     A dictionary containing the edges of every possible path, as created in find_valid_paths
     :return:               A nested dictionary keyed by path, with values corresponding to the symmetric set
                            difference between the edges of those paths."""

    # Note: If you wish to modify the code to assign different scores for each gene node, modifying this function
    # to provide the actual lists of nodes in the SSD might be a good place to start
    # Alternatively, you could provide one symmetric_set_difference_table table per gene node in this function.

    # This is the Symmetric Set Difference table we will be returning
    symmetric_set_difference_table = {}

    path_sets = {}

    # Since we will need the frozenset of each path multiple times, it makes sense to pre-compute them
    for path in path_edges:
        path_sets[path] = frozenset(path_edges[path])

    # To compute that table, we iterate over each combination of paths, and find the length of the symmetric set
    # difference between the frozensets of their path edges.
    for path_a in path_edges:
        symmetric_set_difference_table[path_a] = {}
        for path_b in path_edges:
            symmetric_set_difference_table[path_a][path_b] = len(
                path_sets[path_a].symmetric_difference(path_sets[path_b]))

    # Note: symmetric_set_difference_table[a][b] should equal symmetric_set_difference_table[b][a]. A good optimization
    #  of this function would take advantage of that property to only have to run ~half the computations

    return symmetric_set_difference_table


def build_lossless_path_symmetric_set_difference_table(path_edges):
    """Works like build_path_symmetric_set_difference_table, but returns an SSD filled with all 0s
    :param path_edges:  A dictionary containing the edges of every possible path, as created in find_valid_paths
    :return:            A nested dictionary keyed by path, with values of 0
    """

    # This is the Symmetric Set Difference table we will be returning
    symmetric_set_difference_table = {}

    # To compute that table, we iterate over each combination of paths, and just enter 0.
    for path_a in path_edges.keys():
        symmetric_set_difference_table[path_a] = {}
        for path_b in path_edges.keys():
            symmetric_set_difference_table[path_a][path_b] = 0

    return symmetric_set_difference_table


def build_exit_dicts(dtl_recon_graph):
    """This method goes through every exit event on the reconciliation graph, and builds up three dictionaries that
    pertain to exit events.
    :param dtl_recon_graph: The DTL reconciliaton graph
    :return:    Three dictionaries:
                0. exit_events_by_mapping, which is keyed by mapping node, and where every entry corresponds to a list
                containing every exit event that is a child of that mapping node
                1. exit_events_by_gene, which is keyed by gene node, and where every entry corresponds to a list
                containing every exit event in the group of that gene node
                2. exit_mappings_by_gene, which is keyed by gene node, and where every entry corresponds to a list
                containing every mapping node in group of that gene node, on the condition that that mapping node
                has an exit event
    """
    # These dicts are explained in the docstring directly above here.
    exit_events_by_mapping = {}
    exit_events_by_gene = {}
    exit_mappings_by_gene = {}

    # Here we iterate over every event node of every mapping node.
    for mapping_node in dtl_recon_graph:
        for event in dtl_recon_graph[mapping_node]:

            # Ignore scores and loss events
            if isinstance(event, tuple) and event[0] != 'L':

                gene_node = mapping_node[0]
                # For each of these three blocks, the if statement initializes the list if necessary, and then adds the
                # correct element into the list
                if gene_node not in exit_events_by_gene:
                    exit_events_by_gene[gene_node] = []
                # Use append to avoid adding the elements of the tuple themselves to the list
                exit_events_by_gene[gene_node].append(event)

                if gene_node not in exit_mappings_by_gene:
                    exit_mappings_by_gene[gene_node] = []
                exit_mappings_by_gene[gene_node] += [mapping_node]

                if mapping_node not in exit_events_by_mapping:
                    exit_events_by_mapping[mapping_node] = []
                exit_events_by_mapping[mapping_node].append(event)

    return exit_events_by_mapping, exit_events_by_gene, exit_mappings_by_gene


def intersect_cost(event):
    return 0

def cost(event):
    return 1


def calculate_ancestral_table(species_tree):
    """
    :param species_tree: a species tree, in vertex format and postorder,
    represented as an OrderedDict (output from reformat_tree)
    :return: A nested dictionary. The first dictionary has vertices in the
    tree as keys and the values are dictionaries. These dictionaries have
    as keys vertices of the tree (again) and values which are strings,
    representing how the first index relates to the second (see below
    for info on what certain strings mean). It creates these dictionaries
    by traversing the tree.
    """

    # Initialize the ancestral table which we will be returning
    ancestral_table = dict()

    # Helper dict to help us determine if two nodes are ancestrally related
    descendants = dict()

    # Get all of the vertices in the tree, which are the keys of the species_tree OrderedDict
    vertices = species_tree.keys()

    # Initialize all entries to incomparable to make following calculations easier
    for pair in list(product(vertices, vertices)):  # Cartesian product of all of the vertices

        # Save variables to match format used in previous discussions of this new algorithm
        A = pair[0]
        B = pair[1]

        # Check if we need to make a dictionary for the first vertex
        if A not in ancestral_table:
            ancestral_table[A] = dict()

        # Set all identical pairs to equal while we're in this loop - check for equality here
        if A == B:
            ancestral_table[A][B] = 'eq'

        else:

            # Set all other pairs to incomparable
            ancestral_table[A][B] = 'in'

    # Now loop over all vertex pairs checking for ancestral connections
    for v in vertices:

        # Save the vertices that represent the children into variables
        child1 = species_tree[v][0]
        child2 = species_tree[v][1]

        # Check for leaf nodes
        if child1 is None and child2 is None:
            descendants[v] = []  # Empty list --> no descendants

        else:

            # The descendants of a node are the direct children and those children's children, and so on
            descendants[v] = descendants[child1] + descendants[child2] + [child1] + [child2]

            # Assign relationship between a node and its descendants, and vice versa
            for descendant in descendants[v]:
                ancestral_table[v][descendant] = 'an'
                ancestral_table[descendant][v] = 'des'

    return ancestral_table
  
def calculate_score_double_exit(enter_table, u, gene_tree, uA, dtl_recon_graph_a, uB, dtl_recon_graph_b):
    """This function computes the score of a 'double exit', where both mapping nodes exit immediately."""
    score_double_exit = float('-inf')

    # Test to see if u is a leaf
    if gene_tree[u] == (None, None):
        #TODO: Why do we need to check for contemporaneous events?
        if uA == uB and ('C', (None, None), (None, None)) in dtl_recon_graph_a[uA]:
            score_double_exit = 0
    else:
        uA_exit_events = filter(lambda event: isinstance(event, tuple) and event[0] not in ('C', 'L'),
                                dtl_recon_graph_a[uA])
        uB_exit_events = filter(lambda event: isinstance(event, tuple) and event[0] not in ('C', 'L'),
                                dtl_recon_graph_b[uB])
        for e1 in uA_exit_events:
            child1 = e1[1][0]
            child2 = e1[2][0]
            # B and C are the species nodes of the two mapping nodes of e1
            B = e1[1][1]
            C = e1[2][1]
            for e2 in uB_exit_events:
                # E and F are the species nodes of the two mapping nodes of e2
                # We need to account for the case that the children of u are in opposite order between the two events
                if child1 == e2[1][0]:
                    E = e2[1][1]
                    F = e2[2][1]
                else:
                    E = e2[2][1]
                    F = e2[1][1]
                # Now, we need to turn the species nodes into the correct mapping nodes
                uB = (child1, B)
                uE = (child1, E)
                uC = (child2, C)
                uF = (child2, F)
                # If the score of this iteration's double exit is better than the old one, then the old one will
                # supersede this one
                score_double_exit = max(score_double_exit,
                                        enter_table[child1][uB][uE] + enter_table[child2][uC][uF] \
                                        + (cost(e1) + cost(e2) if e1 != e2 else intersect_cost(0)))

    return score_double_exit


def compute_incomparable_single_exit(enter_table, u, uA, uA_loss_events, uB, uB_loss_events, score_double_exit):

    scores = [score_double_exit]
    for event in uA_loss_events:
        a_child = event[1][1]
        scores += [enter_table[u][(u, a_child)][uB] + cost(event)]
    for event in uB_loss_events:
        b_child = event[1][1]
        scores += [enter_table[u][uA][(u, b_child)] + cost(event)]
    enter_table[u][uA][uB] = max(scores)


def compute_equal_single_exit(enter_table, u, uA, uA_loss_events, uB, uB_loss_events,
                              score_double_exit, single_exit_table):
    assert uA == uB
    scores = [score_double_exit]
    for a_event in uA_loss_events:
        a_child = a_event[1][1]
        for b_event in uB_loss_events:
            b_child = b_event[1][1]
            scores += [enter_table[u][(u, a_child)][(u, b_child)]]

    for event in uA_loss_events:
        a_child = event[1][1]
        scores += [single_exit_table[u][uB][(u, a_child)] + cost(event)]
    for event in uB_loss_events:
        b_child = event[1][1]
        scores += [single_exit_table[u][uA][(u, b_child)] + cost(event)]
    enter_table[u][uA][uB] = max(scores)


def compute_ancestral_single_exit(is_swapped, enter_table, u, uA, uA_loss_events, uB, uB_loss_events,
                              score_double_exit, single_exit_table):
    scores = [score_double_exit]

    # We check to see if which mapping node is the ancestor is swapped from uA an uB to uB an uA. We can't just
    # swap the arguments in that case unfortunately, because enter_table requires the two arguments be entered in the
    # correct direction.
    if not is_swapped:
        for event in uB_loss_events:
            b_child = event[1][1]
            scores += [single_exit_table[u][uA][(u, b_child)] + cost(event)]
        if not uA in single_exit_table[u]:
            single_exit_table[u][uA] = {}
        single_exit_table[u][uA][uB] = max(scores)

        enter_scores = [single_exit_table[u][uA][uB]]
        for event in uA_loss_events:
            a_child = event[1][1]
            enter_scores += [enter_table[u][(u, a_child)][uB] + cost(event)]
        enter_table[u][uA][uB] = max(enter_scores)
    else:
        for event in uA_loss_events:
            a_child = event[1][1]
            scores += [single_exit_table[u][uB][(u, a_child)] + cost(event)]
        if not uB in single_exit_table[u]:
            single_exit_table[u][uB] = {}
        single_exit_table[u][uB][uA] = max(scores)

        enter_scores = [single_exit_table[u][uB][uA]]
        for event in uB_loss_events:
            b_child = event[1][1]
            enter_scores += [enter_table[u][uA][(u, b_child)] + cost(event)]
        enter_table[u][uA][uB] = max(enter_scores)


# Note species tree is assumed to be in vertex format
def new_diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph_a, dtl_recon_graph_b, debug, zero_loss):

    postorder_gene_nodes = list(gene_tree.keys())
    postorder_species_nodes = list(species_tree.keys())
    postorder_group_a = {}
    postorder_group_b = {}
    for u in gene_tree:
        postorder_group_a[u] = filter(lambda mapping: mapping[0] == u, dtl_recon_graph_a)
        postorder_group_a[u] = sorted(postorder_group_a[u], key=lambda mapping: postorder_species_nodes.index(mapping[1]))
        postorder_group_b[u] = filter(lambda mapping: mapping[0] == u, dtl_recon_graph_b)
        postorder_group_b[u] = sorted(postorder_group_b[u], key=lambda mapping: postorder_species_nodes.index(mapping[1]))
    ancestral_table = calculate_ancestral_table(species_tree)
    single_exit_table = {}
    enter_table = {}
    for u in postorder_gene_nodes:
        enter_table[u] = {}
        single_exit_table[u] = {}
        for uA in postorder_group_a[u]:
            enter_table[u][uA] = {}
            for uB in postorder_group_b[u]:
                score_double_exit = calculate_score_double_exit(enter_table, u, gene_tree, uA, dtl_recon_graph_a, uB, dtl_recon_graph_b)
                #print "{0} - {1}".format(uA, uB)
                ancestry = ancestral_table[uA[1]][uB[1]]

                uA_loss_events = filter(lambda event: isinstance(event, tuple) and event[0] == 'L',
                                        dtl_recon_graph_a[uA])
                uB_loss_events = filter(lambda event: isinstance(event, tuple) and event[0] == 'L',
                                        dtl_recon_graph_b[uB])

                # To compute the proper single exit entry, we must know how the two nodes relate to each other. See the
                # header for a more complete explanation on this.
                if ancestry == 'in':
                    compute_incomparable_single_exit(enter_table, u, uA, uA_loss_events, uB, uB_loss_events,
                                                     score_double_exit)
                elif ancestry == 'eq':
                    compute_equal_single_exit(enter_table, u, uA, uA_loss_events, uB, uB_loss_events,
                                              score_double_exit, single_exit_table)
                elif ancestry == 'des':
                    compute_ancestral_single_exit(True, enter_table, u, uA, uA_loss_events, uB, uB_loss_events,
                                                  score_double_exit, single_exit_table)
                elif ancestry == 'an':
                    compute_ancestral_single_exit(False, enter_table, u, uA, uA_loss_events, uB, uB_loss_events,
                                                  score_double_exit, single_exit_table)
                else:
                    raise ValueError("Invalid ancestry type '{0}', check calculate_ancestral_table().".format(ancestry))
        if debug:
            print_table_nicely(enter_table[u], ", ", "EnterTable({0})".format(u))
    # Now, the diameter of this reconciliation will be the maximum entry on the enter table.
    diameter = 0
    for uA in enter_table[gene_tree_root]:
        for uB in enter_table[gene_tree_root][uA]:
            diameter = max(diameter, enter_table[gene_tree_root][uA][uB])

    return diameter

def event_to_string(event):
    return "{0}:{1}{2} {3}{4}".format(str(event[0]), str(event[1][0]), str(event[1][1]),
                                      str(event[2][0]), str(event[2][1]))


def print_table_nicely(table, deliminator, name="\t", dtype="map"):
    """Takes a table (a 2D dict keyed with tuples) and prints a nicely formatted table. Used for debugging and wall art.
    :param table:       The table we wish to print nicely. It is assumed that the rows and columns are exactly the same,
                        and that both the keys and values will fit within 7 characters (room for one tab space)
                        It is also assumed that the keys are tuples of some kind.
    :param deliminator: What string to put in between the elements of the tuples
    :param name:        What this table should be named (upper left)
    :param dtype:        A string corresponding to the type of data. Valid values are 'event', 'path', and 'map'.
    :return:            Nothing, but prints to the screen a lot.
    """

    print ""
    if len(table) > 30:  # Don't spend too long displaying tables.
        print "Table '{1}' is {0}x{0}, which is bigger than the max size of 30.".format(len(table), name)
        return

    line = "\033[4m{0}\033[1m".format(name)  # Underline top row, bold column headers
    for column in table:
        if dtype == "event":
            line += "\t{0}".format(event_to_string(column))
        elif dtype == "path":
            line += "\t{0}{1}{2}{3}{4}".format(column[0][0], column[0][1], deliminator, column[1][0], column[1][1])
        else:
            line += "\t{0}{1}{2}".format(str(column[0]), deliminator, str(column[1]))
    print line + "\033[0m"

    row_num = 0  # Used to alternate row colors

    for row in table:
        row_num += 1
        line_color = "\033[37m" if row_num % 2 == 0 else "\033[0m"

        line = line_color + "\t\033[4m\033[1m"  # Add bolding and underline to row headers
        if dtype == "event":
            line += "{0}".format(event_to_string(row))
        elif dtype == "path":
            line += "{0}{1}{2}{3}{4}".format(row[0][0], row[0][1], deliminator, row[1][0], row[1][1])
        else:
            line += "{0}{1}{2}".format(str(row[0]), deliminator, str(row[1]))
        line += "\033[0m\t" + line_color  # Remove bolding for entries, then return to line color
        for column in table:
            if row == column:
                line += "\033[33m"  # Highlight diagonals
            line += str(table[column][row]) + "\t"
            if row == column:
                line += line_color
        print line
    print "\033[0m"  # Return to default color


def clean_graph(dtl_recon_graph, gene_tree_root):
    """Cleans up the graph created by DTLReconGraph.py by turning removing scores from events and event lists, and removes any
     loss events on the root gene node.
     :param dtl_recon_graph:    The DTL reconciliation graph that we wish to clean
     :param gene_tree_root:     The root of the gene tree of said graph
     :return:                   Nothing, but modifies dtl_recon_graph"""
    for key in dtl_recon_graph:
        # Get rid of all of the random numbers in the event list
        dtl_recon_graph[key] = filter(lambda e: not isinstance(e, (float, int)), dtl_recon_graph[key])
        # The events in the event dtl_recon_graph are stored as lists which cannot be used as dict keys. Let's fix that.
        # The below code is no longer necessary, as it is being filtered out in DTLReconGraph.py

        # for i in range(0, len(dtl_recon_graph[key])):
            # Get rid of the last value, as it is a number we don't need
            # dtl_recon_graph[key][i] = dtl_recon_graph[key][i][:-1]

        # DTLReconGraph should be filtering the loss events on the root node out, so we don't need to worry about it
        # if key[0] == gene_tree_root:
            # dtl_recon_graph[key] = filter(lambda e: not e[0] == 'L', dtl_recon_graph[key])


def write_to_csv(csv_file, costs, filename, mpr_count, diameter, gene_node_count,
                 DTLReconGraph_time_taken, diameter_time_taken):
    """Takes a large amount of information about a diameter solution and appends it as one row to the provided csv file.
    :param csv_file:                    The csv file to write to
    :param costs:                       A string representing the costs used to calculate the DTL recon graph
    :param filename:                    The name of the file that was reconciled
    :param mpr_count:                   The number of MPRS found
    :param diameter:                    The diameter that was found
    :param gene_node_count:             The total number of nodes in the gene tree
    :param DTLReconGraph_time_taken:    The amount of time that DTLReconGraph took to run
    :param diameter_time_taken:         The amount of time that Diameter took to run
    """
    file_exists = os.path.isfile(csv_file)

    with open(csv_file, 'a') as output_file:
        writer = csv.writer(output_file)

        # Write the headers if we need to.
        if not file_exists:
            writer.writerow(["File Name", "Costs", "MPR Count", "Diameter", "Gene Node Count",
                             "DTLReconGraph Computation Time", "Diameter Computation Time", "Date"])

        writer.writerow([filename, costs, mpr_count, diameter, gene_node_count,
                         DTLReconGraph_time_taken, diameter_time_taken, time.strftime("%c")])


def calculate_diameter_from_file(filename, D, T, L, log=None, debug=False, verbose=True):
    """This function computes the diameter of space of MPRs in a DTL reconciliation problem,
     as measured by the symmetric set distance between the events of the two reconciliations of the pair
      that has the highest such difference.

      :param filename:      The path to the newick tree file to reconcile
      :param D:             The cost for duplication events
      :param T:             The cost for transfer events
      :param L:             The cost for loss events
      :param log:           The csv file to output results to (will create it if it does not exist)
      :param debug:         Whether to print out all of the tables
      :return:              Nothing, but we output results to a csv file."""

    # These statements check to make sure that all arguments were entered correctly.
    assert isinstance(log, (str, unicode)) or log is None
    assert isinstance(filename, (str, unicode))
    assert isinstance(D, (int, float))
    assert isinstance(T, (int, float))
    assert isinstance(L, (int, float))
    assert isinstance(debug, bool)

    # Record the time that DTLReconGraph starts
    start_time = time.clock()

    # Get everything we need from DTLReconGraph
    species_tree, gene_tree, dtl_recon_graph, mpr_count, _ = DTLReconGraph.reconcile(filename, D, T, L)

    # Record the time that this code starts

    # The gene tree needs to be in node format, not edge format, so we find that now.
    # (This also puts the gene_tree into postorder, as an ordered dict)
    gene_tree, gene_tree_root, gene_node_count = reformat_tree(gene_tree, "pTop")

    species_tree, species_tree_root, species_node_count = reformat_tree(species_tree, "hTop")

    # The DTL reconciliation graph as provided by DTLReconGraph has some extraneous numbers. We remove those here.
    clean_graph(dtl_recon_graph, gene_tree_root)

    # And record the amount of time DTLReconGraph + cleaning up the graph took
    DTLReconGraph_time_taken = time.clock() - start_time

    if verbose:
        print "Reconciliation Graph Made in \033[33m\033[1m{0} seconds\033[0m".format(DTLReconGraph_time_taken)

    start_time = time.clock()

    # Now we draw the rest of the owl
    diameter = new_diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph, debug, False)

    # And record how long it took to compute the diameter.
    diameter_time_taken = time.clock() - start_time

    #start_time = time.clock()
    #zl_diameter = new_diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph, debug, True)
    #zl_diameter_time_taken = time.clock()-start_time

    if verbose:
        print "The diameter of the given reconciliation graph is \033[33m\033[1m{0}\033[0m, (or \033[33m\033[1m{1}\033[0m if losses do not affect the diameter)".format(diameter, diameter)#zl_diameter)

    # Timing data is inaccurate in debug mode (print statements take too long), so we only give it to the user in non-
    # debug mode.
    if not debug and verbose:
        print "Diameter found in \033[33m\033[1m{0} seconds\033[0m".format(diameter_time_taken)
        print "Total time: \033[33m\033[1m{0} seconds\033[0m".format(diameter_time_taken + DTLReconGraph_time_taken)

    # Now, we write our results to a csv file.
    if log is not None:
        costs = "D: {0} T: {1} L: {2}".format(D, T, L)
        write_to_csv(log + ".csv", costs, filename, mpr_count, diameter, gene_node_count, DTLReconGraph_time_taken, diameter_time_taken)
        #write_to_csv(log + "_zl.csv", costs, filename, mpr_count, zl_diameter, gene_node_count, DTLReconGraph_time_taken,
         #            zl_diameter_time_taken)
    # And we're done.
    return


def repeatedly_calculate_diameter(file_pattern, start, end, d, t, l, log=None, debug=False, verbose=True):
    """Iterates over a lot of input files and finds the diameter of all of them.
    :param file_pattern: A string contains the name of the files to be used, with the counting number replaced with #'s
    :param start:       Numbered file to start on
    :param end:         Numbered file to end with
    :param d:           Duplication event cost
    :param t:           Transfer event cost
    :param l:           Loss event cost
    :param log:         csv file to log results to
    :param debug:       Whether to print out every DP table made (not recommended)
    :return:
    """
    match = re.match("([^#]*)(#+)([^#]*)", file_pattern)
    if not match:
        print "Filepath '" + file_pattern + "' not understood. Please enter the path to your files, with repeated hash marks" \
                                       "(#'s) in place of sequential numbering."
        return
    fill = len(match.group(2))
    if fill < len(str(end-1)) or fill < len(str(start)):
        print "Starting or ending number is larger than '{1}' supports ({0})!".format((10 ** fill) -1, file_pattern)
        return
    print "Running {4} sequential jobs on files '{3}' with DTL of {0},{1},{2}".format(d, t, l, file_pattern, end - start)
    for i in range(start, end):
        cur_file = "{0}{1}{2}".format(match.group(1), str(i).zfill(fill), match.group(3))
        print "Reconciling {0}".format(cur_file)
        try:
            calculate_diameter_from_file(cur_file, d, t, l, log, debug, verbose)
        except IOError:
            print "(File Not Found)"
        except (KeyboardInterrupt, SystemExit):
            raise  # Don't prevent the user from exiting the program.
        except:
            if verbose:
                print traceback.print_exc(sys.exc_traceback)
            print "Could not reconcile file '{0}'. Continuing, but please make sure the file was formatted correctly!"\
                .format(cur_file)


def main():
    """Processes command line arguments"""
    usage = "usage: %prog [options] file d t l"
    p = optparse.OptionParser(usage=usage)
    p.add_option("-l", "--log", dest="logfile", help="writes a logfile in CSV format to LOGFILE", metavar="LOGFILE")
    p.add_option("-i", "--iterate", dest="count", action="store", nargs=2, help="calculates every file matching a "
                                                                                "pattern (defined in the file argument)"
                                                                                " with number MIN to MAX.",
                 metavar="MIN MAX", type=int)
    p.add_option("-d", "--debug", dest="debug", action="store_true", default=False, help= "prints out every DP table with size less"
                                                                           "than 30x30")
    p.add_option("-q", "--quiet", dest="verbose", action="store_false", default=True,
                 help="suppresses (most) text output")

    (options, args) = p.parse_args()
    if len(args) != 4:
        p.error("4 arguments must be provided: file, d, t, and l")
    file = args[0]
    d = float(args[1])
    t = float(args[2])
    l = float(args[3])
    log = options.logfile
    debug = options.debug
    verbose = options.verbose
    if not (log or debug or verbose):
        p.error("some form of output must be specified! (-l or -d must be used when -q is used)")
    elif options.count is not None:
        rep = options.count
        repeatedly_calculate_diameter(file, rep[0], rep[1], d, t, l, log, debug, verbose)
    else:
        calculate_diameter_from_file(file, d, t, l, log, debug, verbose)

if __name__ == "__main__":
    main()

def t(file_name):
    calculate_diameter_from_file(file_name, 2, 3, 1)
