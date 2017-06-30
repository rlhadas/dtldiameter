# Diameter.py
# Written by Eli Zupke and Andrew Ramirez June 2017
# It is based off of Jordan Haack's work on a polynomial-time algorithm for the DTL MPR Diameter Problem

# 0. ON THE PAPER
#
#       This program is based off of an algorithm written by Jordan Haack. There is a paper that describes this
#   algorithm more thoroughly than can be expressed in source code comments. The paper title is "Computing the Diameter
#   of the Space of Maximum Parsimony Reconciliations in the Duplication-Transfer-Loss Model". #TODO add more info
#
#   Note that in this documentation and in that paper, there is a notion of the "group" of a node

# 1. ON TREE REPRESENTATION FORMATS:
#
#       This file deals with trees in two formats: Edge-based formats, and vertex-based formats. The edged-based
#   trees are what are output from DP.py (which returns them straight from newickFormatReader), and are what this
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
#       This file calls the two trees the gene tree and the species tree. Other programs, like DP.py, have different
#   naming conventions (such as "host" for the species tree and "parasite" for the gene tree) because they were coded
#   under different assumptions as to what the two trees represent. Let it be understood that these names are
#   synonymous, and that references to "hTop" or "pTop" refer to the name of the handle of the species and
#   gene trees that DP outputs.


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
#       The tables are described in detail in the technical report, but a summary of their functions follows:
#
#       exit_event_table contains the largest number of event nodes that each pair of reconciliations rooted at events
#   x and y can differ for, where x and y are exit-event nodes in Group(u).
#       exit_mapping_table contains the largest number of event nodes that each pair of reconciliation subtrees rooted
#   at mapping nodes x and y can differ for, on the condition that x and y go immediately to an exit-event.
#       enter_mapping_table is contains the largest number of event nodes that each pair of reconciliation subtrees
#   rooted at x and y can differ for when x and y are used to enter Group(u).
#
#       To get a better sense of how the tables are structured, you can try running an example in debug mode, which
#   will actually print out these tables.

# 5. ON EXIT_X_BY_Y DICTIONARIES:
#
#       There are several helper dictionaries which pertain to exit events. Their names should be pretty self
#   explanatory, but they are described more fully here.
#
#   exit_events_by_mapping, is keyed by mapping node, and every entry corresponds to a list containing every exit event
# that is a child of that mapping node
#   exit_events_by_gene, is keyed by gene node, and where entry corresponds to a list containing every exit event in
# the group of that gene node
#   exit_mappings_by_gene, is keyed by gene node, and where  entry corresponds to a list containing every mapping node
# in group of that gene node, on the condition that that mapping node has an exit event

# -1. DATA STRUCTURE QUICK REFERENCE:
#
#   DTL Reconciliation graph:
#       { mapping_node: [event1, event2, ... eventn, number] ...}
#
#   Mapping node:
#       ('gene_node','SPECIES_NODE')
#   or in loss or contemporary event nodes:
#       (None, None)
#
#   Event node:
#       ('type', child_mapping_node1, child_mapping_node2, number)
#
#   (edge) tree:
#       {('R','N'): ('R','N', ('N','C1'), ('N','C2')) ...}
#       aka:
#       {root_edge: (root_edge[0], root_edge[1], child1_edge, child2_edge) ...}
#
#   vertex_tree:
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
import sys
from collections import OrderedDict


def reformat_tree(tree, root):
    """A recursive function that changes the format of a (species or gene) tree from edge to vertex, as described
    aboce. It returns the tree (in postorder), the root of the tree, and the number of nodes in the tree. The base
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

    next_edges = previous_edges + [current_edge]  # This list contains every node we visited to get here from the current_edge
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


def compute_trivial_exit_event_table(u, exit_event_table):
    """This function computes and stores the score of the exit event on a leaf node 'u' of the gene tree.
    As this event will always be a C event that is shared by all nodes, this value will always be 0.
    :param u:                   The gene node for which we want to compute the table.
    :param exit_event_table:    The dp table we will compute """
    exit_event_table[u] = {}
    exit_event_table[u][('C', (None, None), (None, None))] = {}
    exit_event_table[u][('C', (None, None), (None, None))][('C', (None, None), (None, None))] = 0


def compute_exit_event_table(u, exit_event_table, enter_mapping_table, exit_events_by_gene):
    """This function computes and stores the score of the exit event on a non-leaf node 'u' of the gene tree.
    :param u:                       The gene node for which we want to compute the table.
    :param exit_event_table:        The dp table we will compute
    :param enter_mapping_table:     A dp table that we will use
    :param exit_events_by_gene:
    """

    # Initialize u's exit_event_table entry
    exit_event_table[u] = {}

    # Here we must iterate over every pair of events
    for e1 in exit_events_by_gene[u]:
        child1 = e1[1][0]
        child2 = e1[2][0]
        # B and C are the species nodes of the two mapping nodes of e1
        B = e1[1][1]
        C = e1[2][1]
        exit_event_table[u][e1] = {}
        for e2 in exit_events_by_gene[u]:
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
            exit_event_table[u][e1][e2] = enter_mapping_table[child1][uB][uE] \
                                           + enter_mapping_table[child2][uC][uF] \
                                           + (2 if e1 != e2 else 0)


def compute_exit_mapping_table(u, exit_mapping_table, exit_event_table, exit_mappings_by_gene, exit_events_by_mapping):
    """This function computes and stores the maximum possible score of the exit from gene node u
    :param u:
    :param exit_mapping_table:
    :param exit_event_table:
    :param exit_mappings_by_gene:
    :param exit_events_by_mapping:
    """

    # u_mapping_nodes contains all nodes in Group(u) that have exit events
    u_mapping_nodes = exit_mappings_by_gene[u]

    # Initialize u's exit_mapping_table entry
    exit_mapping_table[u] = {}

    for uA in u_mapping_nodes:
        exit_mapping_table[u][uA] = {}
        for uB in u_mapping_nodes:
            max_value = 0
            for E1 in exit_events_by_mapping[uA]:
                for E2 in exit_events_by_mapping[uB]:
                    max_value = max(exit_event_table[u][E1][E2], max_value)
            exit_mapping_table[u][uA][uB] = max_value


def build_loss_reachable(graph, root_mapping_node):
    """A recursive function to create the list of the mapping nodes that can be reached through loss events.
    (The base case is when you can't get to a loss node from a mapping node)
    :param graph:
    :param root_mapping_node:
    :return: """
    loss_reachable = []

    loss_events = filter(lambda e: e[0] == 'L', graph[root_mapping_node])  # Grab all loss events
    # For explanation of this, see the data structure quick reference.

    for event in loss_events:
        loss_reachable += build_loss_reachable(graph, event[1])
    return loss_reachable + [root_mapping_node]


def compute_enter_mapping_table(u, enter_mapping_table, exit_mapping_table, mapping_node_list, graph, ssd):
    """This function computes the maximum possible score of each pair of mapping nodes for gene node u, and stores each
    one into the enter_mapping_table table for u.
    :param u:
    :param enter_mapping_table:
    :param exit_mapping_table:
    :param mapping_node_list:
    :param graph:
    :param ssd: """

    # TODO: look over this, this is group u, think of removing from fucntion
    u_mapping_nodes = []  # Make a new list that has only the mapping nodes that contain u
    for node in mapping_node_list:
        if node[0] == u:
            u_mapping_nodes.append(node)

    loss_reachable = {}  # Every mapping node (in u)'s list of mapping nodes that can be reached through loss events

    for mapping_node in u_mapping_nodes:
        loss_reachable[mapping_node] = build_loss_reachable(graph, mapping_node)

    enter_mapping_table[u] = {}
    for uA in u_mapping_nodes:
        enter_mapping_table[u][uA] = {}
        for uB in u_mapping_nodes:
            max_score = 0

            # temp_debug_scores = {}  # this table displays the values considered for max_score
            # for map1 in frozenset(loss_reachable[uA] + loss_reachable[uB]):
            #    temp_debug_scores[map1] = {}
            #    for map2 in frozenset(loss_reachable[uA] + loss_reachable[uB]):
            #        temp_debug_scores[map1][map2] = "N/A"

            for uC in loss_reachable[uA]:
                #  Sometimes, we consider values for uC and uD that do not have entries in exit_mapping_table[u].
                #  This means they do not have exit events, and we can ignore them.
                if not uC in exit_mapping_table[u]:
                    break

                for uD in loss_reachable[uB]:

                    if not uD in exit_mapping_table[u][uC]:
                        break
                    score_loss = ssd[(uA[1], uC[1])][(uB[1], uD[1])]
                    score_rest = exit_mapping_table[u][uC][uD]
                    # temp_debug_scores[uC][uD] = score_rest + score_loss
                    max_score = max(max_score, (score_loss + score_rest))

            # print_table_nicely(temp_debug_scores,"","{4}:tmp {0}{1}:{2}{3}".format(uA[0], uA[1], uB[0], uB[1],u))

            enter_mapping_table[u][uA][uB] = max_score


def event_to_string(event):
    return "{0}:{1}{2} {3}{4}".format(str(event[0]), str(event[1][0]), str(event[1][1]),
                                      str(event[2][0]), str(event[2][1]))


def print_table_nicely(table, deliminator, name="\t", type="map"):
    """Takes a table (a 2D dict keyed with tuples) and prints a nicely formatted table. Used for debugging and wall art.
    :param table:
    :param deliminator:
    :param name:
    :param type:
    :return:
    """

    print ""
    if len(table) > 30:  # Don't spend too long displaying tables.
        print "Table '{1}' is {0}x{0}, which is bigger than the max size of 30.".format(len(table), name)
        return

    line = "\033[4m{0}\033[1m".format(name)  # Underline top row, bold column headers
    for column in table:
        if type == "event":
            line += "\t{0}".format(event_to_string(column))
        elif type == "path":
            line += "\t{0}{1}{2}{3}{4}".format(column[0][0], column[0][1], deliminator, column[1][0], column[1][1])
        else:
            line += "\t{0}{1}{2}".format(str(column[0]), deliminator, str(column[1]))
    print line + "\033[0m"

    row_num = 0  # Used to alternate row colors

    for row in table:
        row_num += 1
        line_color = "\033[37m" if row_num % 2 == 0 else "\033[0m"

        line = line_color + "\t\033[4m\033[1m"  # Add bolding and underline to row headers
        if type == "event":
            line += "{0}".format(event_to_string(row))
        elif type == "path":
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
    """Cleans up the graph created by DP.py by turning removing scores from events and event lists, and removes any
     loss events on the root gene node.
     :param dtl_recon_graph:
     :param gene_tree_root: """
    for key in dtl_recon_graph:
        # Get rid of all of the random numbers in the event list
        dtl_recon_graph[key] = filter(lambda e: not isinstance(e, (float, int)), dtl_recon_graph[key])
        # The events in the event dtl_recon_graph are stored as lists which cannot be used as dict keys. Let's fix that.
        #for i in range(0, len(dtl_recon_graph[key])):
            # Get rid of the last value, as it is a number we don't need
            #dtl_recon_graph[key][i] = dtl_recon_graph[key][i][:-1]

        # DP should be filtering the loss events on the root node out, so we don't need to worry about it
        # if key[0] == gene_tree_root:
            # dtl_recon_graph[key] = filter(lambda e: not e[0] == 'L', dtl_recon_graph[key])



def diameter_algorithm(species_tree, vertex_gene_tree, gene_tree_root, dtl_recon_graph, debug, zero_loss):
    """This function is the one that actually computes the diameter. It initializes many dictionaries, computes the
    symmetric set difference between every pair of paths on the species tree, and then runs the algorithm as described
    in Jordan's paper.
    :param species_tree:
    :param vertex_gene_tree:
    :param gene_tree_root:
    :param dtl_recon_graph:
    :param debug:
    :param zero_loss:
    :return: """

    # TODO: Fill this in (all functions)

    # We need to get path_edges, a dict containing all of the edges for each path
    #   and non_trivial_path_list, a subset of path_list with all paths A->A removed.
    path_edges, non_trivial_path_list = find_valid_paths("hTop", species_tree)

    # The dynamic programming tables, as described in section 4 of the header.
    exit_event_table = {}
    exit_mapping_table = {}
    enter_mapping_table = {}

    # This list contains each mapping node in the reconciliation dtl_recon_graph
    all_mapping_nodes = dtl_recon_graph.keys()

    # path_symmetric_set_difference is a 2D dict containing the SSD count for each pair of species node paths.
    if zero_loss:  # If losses don't count for Diameter, then we use an all 0 pSSD table.
        path_symmetric_set_difference = build_lossless_path_symmetric_set_difference_table(path_edges)
    else:
        path_symmetric_set_difference = build_path_symmetric_set_difference_table(path_edges)

    # These dictionaries are explained in section 5 of the header.
    exit_events_by_mapping, exit_events_by_gene, exit_mappings_by_gene = build_exit_dicts(dtl_recon_graph)

    if debug:
        print_table_nicely(path_symmetric_set_difference, "->", "[[SSD]]:")

    # Gene tree's should already be in postorder, as vertex_gene_tree is an ordered dict.
    postorder_gene_vertices = vertex_gene_tree.keys()

    # This next loop is the algorithm just as described in the paper.
    for u in postorder_gene_vertices:

        # Check to see if this is a leaf node of the gene tree
        if vertex_gene_tree[u][0] is not None:

            # If not, then the enter_mapping_table for our children had better exist.
            assert enter_mapping_table[vertex_gene_tree[u][0]] != {} and \
                   enter_mapping_table[vertex_gene_tree[u][1]] != {}

            # Now we must build the exit event table
            compute_exit_event_table(u, exit_event_table, enter_mapping_table, exit_events_by_gene)
        else:

            # If so, u's exit_event_table table is trivial.
            compute_trivial_exit_event_table(u, exit_event_table)

        # Next, we must compute the exit and enter mapping score tables.
        compute_exit_mapping_table(u, exit_mapping_table, exit_event_table, exit_mappings_by_gene, exit_events_by_mapping)
        compute_enter_mapping_table(u, enter_mapping_table, exit_mapping_table, all_mapping_nodes, dtl_recon_graph,
                                    path_symmetric_set_difference)

        # If we are in debug mode, we will want some pretty tables to distract us from the harshness of reality. In
        # real situations, printing to the screen takes too long and the tables are too big to be useful.
        if debug:
            print_table_nicely(exit_event_table[u], ", ", "ExitEventS({0})".format(u), "event")
            print_table_nicely(exit_mapping_table[u], "", "ExitMapS({0})".format(u))
            print_table_nicely(enter_mapping_table[u], "", "EnterMapS({0})".format(u))

    # Now we find the diameter, which will be the total maximum of the gene root's enter table.
    diameter = 0
    for uA in enter_mapping_table[gene_tree_root]:
        for uB in enter_mapping_table[gene_tree_root][uA]:
            diameter = max(enter_mapping_table[gene_tree_root][uA][uB], diameter)
    return diameter


def write_to_csv(csv_file, costs, filename, mpr_count, diameter, gene_node_count,
                 DP_time_taken, diameter_time_taken):
    file_exists = os.path.isfile(csv_file)

    with open(csv_file, 'a') as output_file:
        writer = csv.writer(output_file)

        # Write the headers if we need to.
        if not file_exists:
            writer.writerow(["File Name", "Costs", "MPR Count", "Diameter", "Gene Node Count",
                             "DP Computation Time", "Diameter Computation Time", "Date"])

        writer.writerow([filename, costs, mpr_count, diameter, gene_node_count,
                         DP_time_taken, diameter_time_taken, time.strftime("%c")])


def calculate_diameter_from_file(filename, D, T, L, csv_file="TestLog", debug=False):
    """This function computes the diameter of space of MPRs in a DTL reconciliation problem,
     as measured by the symmetric set distance between the events of the two reconciliations of the pair
      that has the highest such difference.

      :param filename:      The path to the newick tree file to reconcile
      :param D:             The cost for duplication events
      :param T:             The cost for transfer events
      :param L:             The cost for loss events
      :param csv_file:      The csv file to output results to (will create it if it does not exist)
      :param debug:         Whether to print out all of the tables
      :return: Nothing"""

    # These statements check to make sure that all arguments were entered correctly.
    assert isinstance(csv_file, (str, unicode))
    assert isinstance(filename, (str, unicode))
    assert isinstance(D, int)
    assert isinstance(T, int)
    assert isinstance(L, int)
    assert isinstance(debug, bool)

    # Record the time that DP starts
    start_time = time.clock()

    # Get everything we need from DP
    species_tree, gene_tree, dtl_recon_graph, mpr_count = DTLReconGraph.reconcile(filename, D, T, L)

    # And record the amount of time DP took
    DP_time_taken = time.clock() - start_time

    print "Reconciliation Complete in \033[33m\033[1m{0} seconds\033[0m".format(DP_time_taken)

    # Record the time that this code starts
    start_time = time.clock()

    # The gene tree needs to be in node format, not edge format, so we find that now.
    # (This also puts the gene_tree into postorder, as an ordered dict)
    gene_tree, gene_tree_root, gene_node_count = reformat_tree(gene_tree, "pTop")

    # The DTL reconciliation graph as provided by DP has some extraneous numbers. We remove those here.
    clean_graph(dtl_recon_graph, gene_tree_root)

    # Now we draw the rest of the owl
    diameter = diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, debug, False)

    # And record how long it took to compute the diameter.
    diameter_time_taken = time.clock() - start_time

    start_time = time.clock()
    zl_diameter = diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, debug, True)
    zl_diameter_time_taken = time.clock()-start_time


    print "The diameter of the given reconciliation dtl_recon_graph is \033[33m\033[1m{0}\033[0m".format(diameter)

    # Timing data is inaccurate in debug mode (print statements take too long), so we only give it to the user in non-
    # debug mode.
    if not debug:
        print "Diameter found in \033[33m\033[1m{0} seconds\033[0m".format(diameter_time_taken)
        print "Total time: \033[33m\033[1m{0} seconds\033[0m".format(diameter_time_taken + DP_time_taken)

    # Now, we write our results to a csv file.
    costs = "D: {0} T: {1} L: {2}".format(D, T, L)
    write_to_csv(csv_file + ".csv",costs,filename,mpr_count,diameter,gene_node_count,DP_time_taken,diameter_time_taken)
    write_to_csv(csv_file + "_zl.csv", costs, filename, mpr_count, zl_diameter, gene_node_count, DP_time_taken,
                 zl_diameter_time_taken)
    # And we're done.
    return

# -2. COMMAND LINE FUNCTIONS

def rep_calc():
    """Command line function to repeatedly run through numbered files located at TreeLifeData/COG####.newick"""
    if not 7 <= len(sys.argv) <= 8:
        print_help()
        return
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    d = int(sys.argv[4])
    t = int(sys.argv[5])
    l = int(sys.argv[6])
    if len(sys.argv) == 8:
        log = sys.argv[7]
    else:
        log = "Log_File"
    print "Running " + str(end - start) + " sequential jobs on TreeLifeData dataset with DTL of {0},{1},{2}".format(d,t,l)
    for i in range(start, end):
        print "Reconciling COG" + str(i).zfill(4)
        try:
            calculate_diameter_from_file("TreeLifeData/COG{0}.newick".format(str(i).zfill(4)), d, t, l, log, False)
        except IOError:
            print "(File Not Found)"

def print_help():
    """Prints a usage string."""
    print "Usage:"
    print "\ttest: runs a test function"
    print "\tcalc file d t l [logfile]: calculates the diameter of a provided newick file"
    print "\trep start end d t l [logfile]: repeatedly runs calc over the numbered COG files located in TreeLifeData"

def test():
    """Command line function to run a short test."""
    calculate_diameter_from_file("example", 2, 3, 1, "TestLog", True)

def calc():
    """Command line function to calculate the diameter of a file"""
    if not 6 <= len(sys.argv) <= 7:
        print_help()
        return
    file = sys.argv[2]
    d = int(sys.argv[3])
    t = int(sys.argv[4])
    l = int(sys.argv[5])
    if len(sys.argv) == 7:
        log = sys.argv[6]
    else:
        log = "Log_File"
    print "Reconciling " + file
    try:
        calculate_diameter_from_file(file, d, t, l, log, False)
    except IOError:
        print "(File Not Found)"

if __name__ == "__main__":
    """Processes command line arguments"""
    if len(sys.argv) <= 1:
        print_help()
    elif sys.argv[1] in ["-h", "-H", "--help", "--Help"]:
        print_help()
    elif sys.argv[1] == "test":
        test()
    elif sys.argv[1] == "calc":
        calc()
    elif sys.argv[1] == "rep":
        rep_calc()
    else:
        print_help()
