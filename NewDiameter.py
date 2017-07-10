# NewDiameter.py
# Written by Eli Zupke and Andrew Ramirez July 2017
# It is based off of Jordan Haack's more efficient work on a polynomial-time algorithm for the DTL MPR Diameter Problem

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
#       There are two formats in this code-base for trees: Edge-based trees, and vertex-based trees. The edged-based
#   trees are what are output from DTLReconGraph.py (which returns them straight from newickFormatReader). For
#   readability and convenience purposes, both the gene and species trees are converted into vertex-based
#
#   Both formats use dictionaries to store the tree structure, but they differ somewhat
# Example:
#   When N is a node descended from R, with children C1 and C2,
#
#   An entry in the edge-based representation looks like this:
#       {('R','N'): ('R','N', ('N','C1'), ('N','C2')) ...}
#
#   And an entry in the vertex-based representation looks like this:
#       {'N':('C1','C2') ...}
#
#   As this file almost exclusively uses vertex-based trees, vertex-based trees will just be named trees, while edge-
#   based trees contain the word 'edge' in their names.

# 2. ON THE NAMING CONVENTION OF THE TWO TREES:
#
#       This file calls the two trees the gene tree and the species tree. Other programs, like DTLReconGraph.py, have
#   different naming conventions (such as "host" for the species tree and "parasite" for the gene tree) because they
#   were coded under different assumptions as to what the two trees represent. Let it be understood that these names are
#   synonymous, and that references to "hTop" or "pTop" refer to the name of the handle of the species and
#   gene trees that DTLReconGraph outputs.


# 3. ON GENE, SPECIES, AND MAPPING NODES:
#
#       The convention for variable names in this file is to use lowercase letters to represent gene nodes (commonly
#   'u'), uppercase letters to represent species nodes (commonly 'A' and 'B') and pairs of lower and upper case letter
#   to represent mapping nodes (commonly 'uA' and 'uB').


# 4. ON THE DYNAMIC PROGRAMMING TABLES AND THEIR FUNCTIONS:
#
#       The Diameter algorithm involves the use of two dynamic programming tables:
#   The enter_table and the exit_table.
#
#       The tables are described in detail in the paper, but a summary of their functions follows:
#
#       The enter_table has the format [u][uA][uB], and it contains the largest number of event nodes that each pair
#   of reconciliation subtrees rooted at uA and uB can differ for when uA and uB are used to enter Group(u). Running the
#   program in debug mode will print out this table at every u.
#
#       The exit_table has the format [u][uA][uB] (where uA is an ancestor of uB or uA == uB). It contains the  largest
#   number of event nodes that each pair of reconciliation subtrees rooted at uA and uB can differ for when uA leads
#   to an exit event.


# -1. DATA STRUCTURE QUICK REFERENCE:
#
#   Ancestor Comparability:
#       A ancestor of B:    'an'
#       A descendant of B:  'des'
#       A incomparable to B:'in'
#       A is equal to B:    'eq'
#
#   DTL Reconciliation graph (post clean_graph):
#       { mapping_node: [event1, event2, ...] ...}
#
#   Event node (post clean_graph):
#       ('type', child_mapping_node1, child_mapping_node2)
#
#   Mapping node:
#       ('gene_node','SPECIES_NODE')
#   or in loss or contemporary event nodes:
#       (None, None)
#
#   (edge) trees:
#       {('R','N'): ('R','N', ('N','C1'), ('N','C2')) ...}
#       aka:
#       {root_edge: (root_edge[0], root_edge[1], child1_edge, child2_edge) ...}
#
#   vertex_trees:
#       {'N':('C1','C2') ...}

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


def intersect_cost(event):
    """The cost added if both reconciliations being looked at share a particular event"""
    return 0


def cost(event, zero_loss):
    """The cost added if exactly one of the reconciliations being looked at share a particular event."""
    if zero_loss and event[0] == 'L':
        return 0
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

    # Get all of the vertices in the tree
    vertices = [vertex for vertex in species_tree]

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


def calculate_score_both_exit(zero_loss, enter_table, u, gene_tree, uA, dtl_recon_graph_a, uB, dtl_recon_graph_b):
    """This function computes the score of a 'double exit', where both mapping nodes exit immediately."""
    score_both_exit = float('-inf')

    # Test to see if u is a leaf
    if gene_tree[u] == (None, None):
        #TODO: Why do we need to check for contemporaneous events?
        if uA == uB and ('C', (None, None), (None, None)) in dtl_recon_graph_a[uA]:
            score_both_exit = 0
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

                score_both_exit = max(score_both_exit,
                                        enter_table[child1][uB][uE] + enter_table[child2][uC][uF] +
                                        (cost(e1, zero_loss) + cost(e2, zero_loss) if e1 != e2
                                           else intersect_cost(0)))

    return score_both_exit


def calculate_incomparable_enter_score(zero_loss, enter_table, u, uA, uA_loss_events, uB, uB_loss_events,
                                       score_both_exit):
    """Returns the enter table entry for [uA][uB] with the assumption that A is on a different part of the species
    tree from B
    :param zero_loss:           Whether losses should not count
    :param enter_table:         The DP table we are computing part of
    :param u:                   The gene node whose group we are in
    :param uA:                  The first mapping node to compare
    :param uA_loss_events:      A list of the loss events on that mapping node
    :param uB:                  The second mapping node to compare
    :param uB_loss_events:      A list of the loss events on that mapping node
    :param score_both_exit:   The score of the double-exit that was previously calculated for uA and uB
    """
    scores = [score_both_exit]
    for event in uA_loss_events:
        a_child = event[1][1]
        scores += [enter_table[u][(u, a_child)][uB] + cost(event, zero_loss)]
    for event in uB_loss_events:
        b_child = event[1][1]
        scores += [enter_table[u][uA][(u, b_child)] + cost(event, zero_loss)]
    return max(scores)


def calculate_equal_enter_score(zero_loss, enter_table, u, uA, uA_loss_events, uB, uB_loss_events,
                                score_both_exit, exit_table_a, exit_table_b):
    """Returns the enter table entry for [uA][uB] with the assumption that uA equals uB (but they might have different
    loss events leaading from them!
    :param zero_loss:           Whether losses should not count
    :param enter_table:         The DP table we are computing part of
    :param u:                   The gene node whose group we are in
    :param uA:                  The first mapping node to compare
    :param uA_loss_events:      A list of the loss events on that mapping node
    :param uB:                  The second mapping node to compare
    :param uB_loss_events:      A list of the loss events on that mapping node
    :param score_both_exit:   The score of the double-exit that was previously calculated for uA and uB
    :param exit_table:   The single exit table, which contains information about the single exit events for
                                the mapping nodes' children.
    """
    # If uA does not equal uB, then something's gone horribly wrong.
    assert uA == uB

    # Build up a list of the possible scores of this pair of mapping nodes, so that we can find the maximum later.
    scores = [score_both_exit]
    for a_event in uA_loss_events:
        a_child = a_event[1][1]
        for b_event in uB_loss_events:
            b_child = b_event[1][1]
            scores += [enter_table[u][(u, a_child)][(u, b_child)]]

    for event in uA_loss_events:
        a_child = event[1][1]
        scores += [exit_table_b[u][uB][(u, a_child)] + cost(event, zero_loss)]
    for event in uB_loss_events:
        b_child = event[1][1]
        scores += [exit_table_a[u][uA][(u, b_child)] + cost(event, zero_loss)]
    return max(scores)


def calculate_ancestral_enter_score(zero_loss, is_swapped, enter_table, u, uA, uA_loss_events, uB, uB_loss_events,
                                    score_both_exit, exit_table_a, exit_table_b):
    """Returns the enter table entry for [uA][uB] with the assumption that A is an ancestor of B (if is_swapped is
    false) or that B is an ancestor of A (if is_swapped is true). In both cases, it will compute the single exit
    table entry of the pair (with the ancestor going first, of course).
    :param zero_loss:           Whether losses should not count
    :param is_swapped:          Whether B is an ancestor of A (instead of the assumed A is an ancestor of B)
    :param enter_table:         The DP table we are computing part of
    :param u:                   The gene node whose group we are in
    :param uA:                  The first mapping node to compare
    :param uA_loss_events:      A list of the loss events on that mapping node
    :param uB:                  The second mapping node to compare
    :param uB_loss_events:      A list of the loss events on that mapping node
    :param score_both_exit:   The score of the double-exit that was previously calculated for uA and uB
    :param exit_table:   The single exit table, which contains information about the single exit events for
                                the mapping nodes' children."""

    # In both cases, we will need to tally up the scores of any loss events on the descendant. Scores will hold those
    # values, and the score of a double exit.
    scores = [score_both_exit]

    # We check to see if which mapping node is the ancestor is swapped from uA an uB to uB an uA. We can't just
    # swap the arguments in that case unfortunately, because enter_table requires the two arguments be entered in the
    # correct direction.
    if not is_swapped:
        # uA is an ancestor to uB
        # Tally up the scores of the descendant's (uB's) loss events
        for event in uB_loss_events:
            b_child = event[1][1]
            # Add the score of taking this loss (the exit_table's entry for the mapping node that this loss
            # leads to, plus the cost of a loss)
            scores += [exit_table_a[u][uA][(u, b_child)] + cost(event, zero_loss)]

        # Initialize the ancestor's (uA) entry in exit_table, if need be.
        if uA not in exit_table_a[u]:
            exit_table_a[u][uA] = {}
        exit_table_a[u][uA][uB] = max(scores)

        enter_scores = [exit_table_a[u][uA][uB]]
        for event in uA_loss_events:
            a_child = event[1][1]
            enter_scores += [enter_table[u][(u, a_child)][uB] + cost(event, zero_loss)]
        return max(enter_scores)
    else:
        # uB is an ancestor to uA
        # Tally up the scores of the descendant's (uA's) loss events
        for event in uA_loss_events:
            a_child = event[1][1]
            # Add the score of taking this loss (the exit_table's entry for the mapping node that this loss
            # leads to, plus the cost of a loss)
            scores += [exit_table_b[u][uB][(u, a_child)] + cost(event, zero_loss)]

        # Initialize the ancestor's (uB) entry in exit_table, if need be.
        if uB not in exit_table_b[u]:
            exit_table_b[u][uB] = {}
        exit_table_b[u][uB][uA] = max(scores)
        
        enter_scores = [exit_table_b[u][uB][uA]]
        for event in uB_loss_events:
            b_child = event[1][1]
            enter_scores += [enter_table[u][uA][(u, b_child)] + cost(event, zero_loss)]
        return max(enter_scores)


def new_diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph_a, dtl_recon_graph_b, debug, zero_loss):
    """
    This function finds the diameter of a reconciliation graph, as measured by the largest symmetric set difference
     of any two reconciliation trees inside of a reconciliation graph. While you can get standard diameter behaviour
     by making dtl_recon_graph_a equal dtl_recon_graph_b, arbitrary restrictions may be placed on which nodes are
     selected by choosing different graphs, for example by limiting one of the graphs to a single reconciliation tree
     to find that tree's distance to the furthest reconciliation.
    :param species_tree:        The species tree (in vertex form)
    :param gene_tree:           The gene tree (in vertex form)
    :param gene_tree_root:      The root of the gene tree
    :param dtl_recon_graph_a:   One of the two DTL reconcilation graphs to make the diameter from.
    :param dtl_recon_graph_b:   The other reconciliation graph. Both must share the same species and gene trees.
    :param debug:               Whether or not to print out pretty tables
    :param zero_loss:           Whether losses should count at all
    :return:                    The diameter of the reconciliation.
    """

    postorder_gene_nodes = list(gene_tree.keys())
    postorder_species_nodes = list(species_tree.keys())
    postorder_group_a = {}
    postorder_group_b = {}
    for u in gene_tree:
        # First we make the dictionary only contain nodes that have this gene node
        postorder_group_a[u] = filter(lambda mapping: mapping[0] == u, dtl_recon_graph_a)
        # Then we sort the dictionary into postorder, using the species node's index in the postorder species list as a
        # guide.
        postorder_group_a[u] = sorted(postorder_group_a[u], key=lambda mapping: postorder_species_nodes.index(mapping[1]))
        # And we do it again for group B
        # First we make the dictionary only contain nodes that have this gene node
        postorder_group_b[u] = filter(lambda mapping: mapping[0] == u, dtl_recon_graph_b)
        # Then we sort the dictionary into postorder, using the species node's index in the postorder species list as a
        # guide.
        postorder_group_b[u] = sorted(postorder_group_b[u], key=lambda mapping: postorder_species_nodes.index(mapping[1]))
    ancestral_table = calculate_ancestral_table(species_tree)

    print_table_nicely(ancestral_table, ", ", "Ancestral", "literal")

    exit_table_a = {}
    exit_table_b = {}

    enter_table = {}

    for u in postorder_gene_nodes:
        enter_table[u] = {}
        exit_table_a[u] = {}
        exit_table_b[u] = {}
        for uA in postorder_group_a[u]:
            enter_table[u][uA] = {}
            for uB in postorder_group_b[u]:
                score_both_exit = calculate_score_both_exit(zero_loss, enter_table, u, gene_tree, uA, dtl_recon_graph_a, uB, dtl_recon_graph_b)


                #print "{0} - {1}".format(uA, uB)
                ancestry = ancestral_table[uA[1]][uB[1]]

                uA_loss_events = filter(lambda event: isinstance(event, tuple) and event[0] == 'L',
                                        dtl_recon_graph_a[uA])
                uB_loss_events = filter(lambda event: isinstance(event, tuple) and event[0] == 'L',
                                        dtl_recon_graph_b[uB])

                # To compute the proper single exit entry, we must know how the two nodes relate to each other. See the
                # header for a more complete explanation on this data structure.
                if ancestry == 'in':
                    score = calculate_incomparable_enter_score(zero_loss, enter_table, u, uA, uA_loss_events, uB, uB_loss_events,
                                                               score_both_exit)
                elif ancestry == 'eq':
                    score = calculate_equal_enter_score(zero_loss, enter_table, u, uA, uA_loss_events, uB, uB_loss_events,
                                                        score_both_exit, exit_table_a, exit_table_b)
                # The only difference between the 'des' and 'an' cases are whether the nodes should be swapped
                elif ancestry == 'des':
                    score = calculate_ancestral_enter_score(zero_loss, True, enter_table, u, uA, uA_loss_events, uB, uB_loss_events,
                                                            score_both_exit, exit_table_a, exit_table_b)
                elif ancestry == 'an':
                    score = calculate_ancestral_enter_score(zero_loss, False, enter_table, u, uA, uA_loss_events, uB, uB_loss_events,
                                                            score_both_exit, exit_table_a, exit_table_b)
                else:
                    raise ValueError("Invalid ancestry type '{0}', check calculate_ancestral_table().".format(ancestry))
                enter_table[u][uA][uB] = score
                if debug:
                    print "{0} -{1}-> {2}, Double-equal\t{3}\tScore:{4}".format(uA,ancestry,uB,score_both_exit,score)


        if debug:
            print_table_nicely(enter_table[u], ", ", "EnterTable({0})".format(u))

    if debug:
        print "Exit Table A: {0}".format(exit_table_a)
        print ""
        print "Exit Table b: {0}".format(exit_table_b)
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
    for column in table[table.keys()[0]]:
        if dtype == "event":
            line += "\t{0}".format(event_to_string(column))
        elif dtype == "literal":
            line+=  "\t{0}".format(column)
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
        elif dtype == "literal":
            line += "{0}".format(column)
        else:
            line += "{0}{1}{2}".format(str(row[0]), deliminator, str(row[1]))
        line += "\033[0m\t" + line_color  # Remove bolding for entries, then return to line color
        for column in table[row]:
            if row == column:
                line += "\033[33m"  # Highlight diagonals
            line += str(table[row][column]) + "\t"
            if row == column:
                line += line_color
        print line
    print "\033[0m"  # Return to default color


def clean_graph(dtl_recon_graph, gene_tree_root):
    """Cleans up the graph created by DTLReconGraph.py by turning removing scores from event lists
     :param dtl_recon_graph:    The DTL reconciliation graph that we wish to clean
     :param gene_tree_root:     The root of the gene tree of said graph
     :return:                   Nothing, but modifies dtl_recon_graph"""
    for key in dtl_recon_graph:
        # Get rid of all of the random numbers in the event list
        dtl_recon_graph[key] = filter(lambda e: not isinstance(e, (float, int)), dtl_recon_graph[key])


def write_to_csv(csv_file, costs, filename, mpr_count, diameter, gene_node_count, species_node_count,
                 DTLReconGraph_time_taken, diameter_time_taken):
    """Takes a large amount of information about a diameter solution and appends it as one row to the provided csv file.
    :param csv_file:                    The csv file to write to
    :param costs:                       A string representing the costs used to calculate the DTL recon graph
    :param filename:                    The name of the file that was reconciled
    :param mpr_count:                   The number of MPRS found
    :param diameter:                    The diameter that was found
    :param gene_node_count:             The total number of nodes in the gene tree
    :param species_node_count:          The total number of nodes in the species tree
    :param DTLReconGraph_time_taken:    The amount of time that DTLReconGraph took to run
    :param diameter_time_taken:         The amount of time that Diameter took to run
    """
    file_exists = os.path.isfile(csv_file)

    with open(csv_file, 'a') as output_file:
        writer = csv.writer(output_file)

        # Write the headers if we need to.
        if not file_exists:
            writer.writerow(["File Name", "Costs", "MPR Count", "Diameter", "Gene Node Count", "Species Node Count",
                             "DTLReconGraph Computation Time", "Diameter Computation Time", "Date"])

        writer.writerow([filename, costs, mpr_count, diameter, gene_node_count, species_node_count,
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
    edge_species_tree, edge_gene_tree, dtl_recon_graph, menpmn, mdenpmn, data, mpr_count, _ = DTLReconGraph.reconcile(filename, D, T, L)

    # Record the time that this code starts

    # The gene tree needs to be in node format, not edge format, so we find that now.
    # (This also puts the gene_tree into postorder, as an ordered dict)
    gene_tree, gene_tree_root, gene_node_count = reformat_tree(edge_gene_tree, "pTop")

    species_tree, species_tree_root, species_node_count = reformat_tree(edge_species_tree, "hTop")

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

    start_time = time.clock()
    zl_diameter = new_diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph, debug, True)
    zl_diameter_time_taken = time.clock()-start_time

    if verbose:
        print "The diameter of the given reconciliation graph is \033[33m\033[1m{0}\033[0m, (or \033[33m\033[1m{1}\033[0m if losses do not affect the diameter)".format(diameter, zl_diameter)

    # Timing data is inaccurate in debug mode (print statements take too long), so we only give it to the user in non-
    # debug mode.
    if not debug and verbose:
        print "Diameter found in \033[33m\033[1m{0} seconds\033[0m".format(diameter_time_taken)
        print "Total time: \033[33m\033[1m{0} seconds\033[0m".format(diameter_time_taken + DTLReconGraph_time_taken)

    # Now, we write our results to a csv file.
    if log is not None:
        costs = "D: {0} T: {1} L: {2}".format(D, T, L)
        write_to_csv(log + ".csv", costs, filename, mpr_count, diameter, gene_node_count, species_node_count, DTLReconGraph_time_taken, diameter_time_taken)
        write_to_csv(log + "_zl.csv", costs, filename, mpr_count, zl_diameter, gene_node_count, species_node_count, DTLReconGraph_time_taken,
                     zl_diameter_time_taken)
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


def t(file_name):
    calculate_diameter_from_file(file_name, 2, 3, 1)


if __name__ == "__main__":
    main()
