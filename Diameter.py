import DP
import time
import csv
import os.path
from collections import OrderedDict

def reformat_tree(tree, root):
    """A recursive function that changes the format of a (species or gene) tree from edge to vertex, for example:
    ('A','B'): ('A','B',('B',C1),('B',C2)) would become 'B':(C1,C2). It returns the tree (in postorder), the root of
    the tree, and the number of nodes in the tree. The base case of this function is when """

    # This line catches the "xTop" handle and replaces
    new_root = root[1] if isinstance(root, tuple) else tree[root][1]

    child1 = tree[root][2][1] if tree[root][2] is not None else None  # These lines handle the leaves, where
    child2 = tree[root][3][1] if tree[root][3] is not None else None  # there is None in the place of a tuple

    # This is the tree that we will be returning. We will add the subtrees of our children first, then add this node.
    new_tree = OrderedDict()  # This has to be an OrderedDict, otherwise we can't guarantee it's in postorder

    # This is the number of nodes in the subtree rooted at each of our children
    child1_count = 0
    child2_count = 0

    if child1 is not None:  # If this node has children, then we need to add their children's subtrees to the dict
        child1_tree, _, child1_count = reformat_tree(tree, tree[root][2])
        new_tree.update(child1_tree)
    if child2 is not None:
        child2_tree, _, child2_count = reformat_tree(tree, tree[root][3])
        new_tree.update(child2_tree)
    new_tree.update({new_root: (child1, child2)})  # We add this node last, to ensure postorderness.

    return new_tree, new_root, (child1_count + child2_count + 1)


def find_valid_paths(root, previous_edges, tree):
    """A recursive algorithm to find all of the valid paths through a binary 'tree' rooted at 'root', as well as a dict
     containing the nodes in each path."""

    next_edges = previous_edges + [root]  # This list contains every node we visited to get here from the root
    paths = []
    path_edges = {}  # This is the edges contained in each path.
    for i in range(0, len(next_edges)):
        source_node = next_edges[i]
        # new_path becomes every path that ends in this value, including A->A
        if i >= len(next_edges)-1:
            new_path = (tree[source_node][1], tree[source_node][1])  # An ugly hack to account for the top of the tree
        else:
            new_path = (tree[source_node][1], root[1])  # This is what we should always be able to do
        paths += [new_path]
        path_edges[new_path] = next_edges[i+1:len(next_edges)]  # This is the list of edges

    child1 = tree[root][2]
    child2 = tree[root][3]

    if child1 is not None:  # Then this Node is not a leaf Node, so we need to add this Node's children
        child1_paths, child1_path_nodes = find_valid_paths(child1, next_edges, tree)
        child2_paths, child2_path_nodes = find_valid_paths(child2, next_edges, tree)
        paths += child1_paths + child2_paths
        child1_path_nodes.update(child2_path_nodes)
        path_edges.update(child1_path_nodes)
        # Otherwise, we have reached the end of the tree (the base case)
    return paths, path_edges


def find_path_lists(root, tree):
    """A function that uses find_valid_paths"""

    paths, path_edges = find_valid_paths(root, [], tree)

    non_trivial_path_list = list(paths)  # Strip the trivial paths (A->A) because they will not have loss events
    for path in non_trivial_path_list:
        if path[0] == path[1]:
            non_trivial_path_list.remove(path)
    return paths, path_edges, non_trivial_path_list


def compute_path_symmetric_set_difference_table(path_list, path_edges):
    """Computes the table containing the number of nodes in the symmetric set difference between any two paths on the
     species tree sTree. This is used in assigning a score to two paths' losses."""

    # Note: If you wish to modify the code to assign different scores for each gene node, modifying this function
    # to provide the actual lists of nodes in the SSD might be a good place to start
    # Alternatively, you could provide one ssd table per gene node in this function.

    ssd = {}  # This is the Symmetric Set Difference table we will be returning

    for path_a in path_list:
        a_nodes = frozenset(path_edges[path_a])
        ssd[path_a] = {}
        for path_b in path_list:
            b_nodes = frozenset(path_edges[path_b])
            ssd[path_a][path_b] = len(a_nodes.symmetric_difference(b_nodes))

    return ssd

def compute_lossles_path_symmetric_set_difference_table(path_list, path_edges):
    """Works like compute_path_symmetric_set_difference_table, but returns an SSD filled with all 0s"""

    ssd = {}  # This is the Symmetric Set Difference table we will be returning

    for path_a in path_list:
        ssd[path_a] = {}
        for path_b in path_list:
            ssd[path_a][path_b] = 0

    return ssd


def build_exit_dicts(graph):
    """Builds a dict containing lists of exit event nodes, each list keyed by the gene node."""
    # The graph returned by DP is keyed by mapping nodes, so to get the gene node we take the first element of the key.
    exit_event_graph = {}  # Every mapping node's list of exit events
    exit_event_dict = {}  # Every gene node's list of exit events
    exit_mapping_dict = {}  # Every gene node's list of mapping nodes that have exit events
    for mapping_node in graph:
        for event in graph[mapping_node]:
            if isinstance(event, tuple) and event[0] != 'L':  # Ignore scores and loss events
                gene_node = mapping_node[0]
                if gene_node not in exit_event_dict:
                    exit_event_dict[gene_node] = []
                exit_event_dict[gene_node].append(event)

                if gene_node not in exit_mapping_dict:
                    exit_mapping_dict[gene_node] = []
                exit_mapping_dict[gene_node] += [mapping_node]

                if mapping_node not in exit_event_graph:
                    exit_event_graph[mapping_node] = []
                exit_event_graph[mapping_node].append(event)

    return exit_event_graph, exit_event_dict, exit_mapping_dict


def compute_trivial_exit_event_table(u, exit_event_scores):
    """This function computes and stores the score of the exit event on a leaf node 'u' of the gene tree.
    As this event will always be a C event that is shared by all nodes, this value will always be 0."""
    exit_event_scores[u] = {}
    exit_event_scores[u][('C', (None, None), (None, None))] = {}
    exit_event_scores[u][('C', (None, None), (None, None))][('C', (None, None), (None, None))] = 0


def compute_exit_event_table(u, exit_event_scores, enter_mapping_scores, exit_event_list):
    """This function computes and stores the score of the exit event on a non-leaf node 'u' of the gene tree."""

    for e1 in exit_event_list[u]:
        child1 = e1[1][0]
        child2 = e1[2][0]
        # B and C are the species nodes of the two mapping nodes of e1
        B = e1[1][1]
        C = e1[2][1]
        exit_event_scores[u][e1] = {}
        for e2 in exit_event_list[u]:
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
            exit_event_scores[u][e1][e2] = enter_mapping_scores[child1][uB][uE] \
                                           + enter_mapping_scores[child2][uC][uF] \
                                           + (2 if e1 != e2 else 0)


def compute_exit_mapping_table(u, exit_mapping_scores, exit_event_scores, exit_mapping_node_dict, exit_event_graph):
    """This function computes and stores the maximum possible score of the exit from gene node u"""

    print_table_nicely(exit_event_scores[u],"-",str(u),"event")

    u_mapping_nodes = exit_mapping_node_dict[u]
    exit_mapping_scores[u] = {}
    for uA in u_mapping_nodes:
        exit_mapping_scores[u][uA] = {}
        for uB in u_mapping_nodes:
            max_value = 0
            for E1 in exit_event_graph[uA]:
                for E2 in exit_event_graph[uB]:
                    max_value = max(exit_event_scores[u][E1][E2], max_value)
            exit_mapping_scores[u][uA][uB] = max_value


def create_loss_reachable(graph, root_mapping_node):
    """A recursive function to create the list of the mapping nodes that can be reached through loss events.
    (The base case is when you can't get to a loss node from a mapping node)"""
    loss_reachable = []

    loss_events = filter(lambda e: e[0] == 'L', graph[root_mapping_node])  # Grab all loss events
    for Event in loss_events:
        loss_reachable += create_loss_reachable(graph, Event[1])
    return loss_reachable + [root_mapping_node]


def compute_enter_mapping_table(u, enter_mapping_scores, exit_mapping_scores, mapping_node_list, graph, ssd):
    """This function computes the maximum possible score of each pair of mapping nodes for gene node u, and stores each
    one into the enter_mapping_scores table for u."""

    u_mapping_nodes = []  # Make a new list that has only the mapping nodes that contain u
    for node in mapping_node_list:
        if node[0] == u:
            u_mapping_nodes.append(node)

    loss_reachable = {}  # Every mapping node (in u)'s list of mapping nodes that can be reached through loss events

    for mapping_node in u_mapping_nodes:
        loss_reachable[mapping_node] = create_loss_reachable(graph, mapping_node)

    enter_mapping_scores[u] = {}
    for uA in u_mapping_nodes:
        enter_mapping_scores[u][uA] = {}
        for uB in u_mapping_nodes:
            max_score = 0

            #temp_debug_scores = {}  # this table displays the values considered for max_score
            #for map1 in frozenset(loss_reachable[uA] + loss_reachable[uB]):
            #    temp_debug_scores[map1] = {}
            #    for map2 in frozenset(loss_reachable[uA] + loss_reachable[uB]):
            #        temp_debug_scores[map1][map2] = "N/A"

            for uC in loss_reachable[uA]:
                #  Sometimes, we consider values for uC and uD that do not have entries in exit_mapping_scores[u].
                #  This means they do not have exit events, and we can ignore them.
                if not uC in exit_mapping_scores[u]:
                    break

                for uD in loss_reachable[uB]:

                    if not uD in exit_mapping_scores[u][uC]:
                        break
                    score_loss = ssd[(uA[1], uC[1])][(uB[1], uD[1])]
                    score_rest = exit_mapping_scores[u][uC][uD]
                    #temp_debug_scores[uC][uD] = score_rest + score_loss
                    max_score = max(max_score, (score_loss + score_rest))

            #print_table_nicely(temp_debug_scores,"","{4}:tmp {0}{1}:{2}{3}".format(uA[0], uA[1], uB[0], uB[1],u))

            enter_mapping_scores[u][uA][uB] = max_score


def event_to_string(event):
    return "{0}:{1}{2} {3}{4}".format(str(event[0]), str(event[1][0]), str(event[1][1]),
                                      str(event[2][0]), str(event[2][1]))


def print_table_nicely(table, deliminator, name="\t", type="map"):
    """Takes a table (a 2D dict keyed with tuples) and prints a nicely formatted table. Used for debugging."""

    print ""
    if len(table) > 30:  # Don't spend too long displaying tables.
        print "Table '{1}' is {0}x{0}, which is bigger than the max size of 30.".format(len(table),name)
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
            line += "{0}{1}{2}{3}{4}".format(row[0][0],row[0][1], deliminator, row[1][0], row[1][1])
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


def clean_graph(graph, gene_tree_root):
    """Cleans up the graph created by DP.py by turning removing scores from events and event lists, and removes any
     loss events on the root gene node."""

    for key in graph:
        # Get rid of all of the random numbers in the event list
        graph[key] = filter(lambda e: not isinstance(e, (float, int)), graph[key])
        # The events in the event graph are stored as lists which cannot be used as dict keys. Let's fix that.
        for i in range(0, len(graph[key])):
            # Get rid of the last value, as it is a number we don't need
            graph[key][i] = graph[key][i][:-1]
        if key[0] == gene_tree_root:
            graph[key] = filter(lambda e: not e[0] == 'L', graph[key])


def diameter_algorithm(edge_species_tree, gene_tree, gene_tree_root, graph, debug, zero_loss):
    """This function is the one that actually computes the diameter. It initializes many dictionaries, computes the
    symmetric set difference between every pair of paths on the species tree, and then runs the algorithm as described
    in Jordan's paper."""

    # We need to get path_list, a list of all of the valid paths in the species tree with format (src, dest)
    #   path_edges, a dict containing all of the edges for each path
    #   and non_trivial_path_list, a subset of path_list with all paths A->A removed.
    path_list, path_edges, non_trivial_path_list = find_path_lists("hTop", edge_species_tree)

    # The key format for the next three dicts are as follows: [u][x][y], where u is a node on the gene tree, and
    # x and y are either event nodes (represented as E1 and E2) or mapping nodes (represented as uA and uB).

    # exit_event_scores is a dict containing the largest number of event nodes that each pair of reconciliation subtrees
    # rooted at E1 and E2 can differ for, where E1 and E2 are exit-event nodes in Group(u).
    exit_event_scores = {}

    # exit_mapping_scores is a dict containing the largest number of event nodes that each pair of reconciliation
    # subtrees rooted at uA and uB can differ for, on the condition that uA and uB go immediately to an exit-event.
    exit_mapping_scores = {}

    # enter_mapping_scores is a dict containing the largest number of event nodes that each pair of reconciliation
    # subtrees rooted at uA and uB can differ for when uA and uB are used to enter Group(u).
    enter_mapping_scores = {}

    # This list contains each mapping node in the reconciliation graph
    mapping_node_list = graph.keys()

    # path_symmetric_set_difference is a 2D dict containing the SSD count for each pair of species node paths.
    if zero_loss:  # If losses don't count for Diameter, then we use an all 0 pSSD table.
        path_symmetric_set_difference = compute_lossles_path_symmetric_set_difference_table(path_list,path_edges)
    else:
        path_symmetric_set_difference = compute_path_symmetric_set_difference_table(path_list, path_edges)

    # exit_event_graph is a dict that contains each exit (non-loss) event in the reconciliation graph, keyed by mapping node.
    #   exit_event_dict is the same as exit_event_graph except it is keyed by gene node.
    #   exit_mapping_node_dict is a dict that contains every mapping node that contains an exit event, keyed by gene node.
    exit_event_graph, exit_event_dict, exit_mapping_node_dict = build_exit_dicts(graph)

    if debug:
        print_table_nicely(path_symmetric_set_difference, "->", "[[SSD]]:")

    # We need to initialize these dicts for every gene node.
    postorder_gene_vertices = gene_tree.keys()
    for u in postorder_gene_vertices:
        enter_mapping_scores[u] = {}
        exit_mapping_scores[u] = {}
        exit_event_scores[u] = {}

    # This next loop is the algorithm just as described in the paper.
    for u in postorder_gene_vertices:

        # Check to see if this is a leaf node of the gene tree
        if gene_tree[u][0] is not None:

            # If not, then the enter_mapping_scores for our children had better exist.
            assert enter_mapping_scores[gene_tree[u][0]] != {} and \
                   enter_mapping_scores[gene_tree[u][1]] != {}

            # Now we must build the exit event table
            compute_exit_event_table(u, exit_event_scores, enter_mapping_scores, exit_event_dict)
        else:

            # If so, u's exit_event_scores table is trivial.
            compute_trivial_exit_event_table(u, exit_event_scores)

        # Next, we must compute the exit and enter mapping score tables.
        compute_exit_mapping_table(u, exit_mapping_scores, exit_event_scores, exit_mapping_node_dict, exit_event_graph)
        compute_enter_mapping_table(u, enter_mapping_scores, exit_mapping_scores, mapping_node_list, graph,
                                    path_symmetric_set_difference)

        # If we are in debug mode, we will want some pretty tables to distract us from the harshness of reality. In
        # real situations, printing to the screen takes too long and the tables are too big to be useful.
        if debug:
            print_table_nicely(exit_event_scores[u], ", ", "ExitEventS({0})".format(u), "event")
            print_table_nicely(exit_mapping_scores[u], "", "ExitMapS({0})".format(u))
            print_table_nicely(enter_mapping_scores[u], "", "EnterMapS({0})".format(u))

    # Now we find the diameter, which will be the total maximum of the gene root's enter table.
    diameter = 0
    for uA in enter_mapping_scores[gene_tree_root]:
        for uB in enter_mapping_scores[gene_tree_root][uA]:
            diameter = max(enter_mapping_scores[gene_tree_root][uA][uB], diameter)
    return diameter


def calculate_diameter_from_file(filename, D, T, L, csv_file="TestLog.csv", debug=False, zero_loss=False):
    """This function computes the diameter of space of MPRs in a DTL reconciliation problem,
     as measured by the symmetric set distance between the events of the two reconciliations of the pair
      that has the highest such difference.

      :param filename:      The path to the newick tree file to reconcile
      :param D:             The cost for duplication events
      :param T:             The cost for transfer events
      :param L:             The cost for loss events
      :param csv_file:      The csv file to output results to (will create it if it does not exist)
      :param debug:         Whether to print out all of the tables
      :param zero_loss:     Whether the diameter should NOT include loss events
      :return: Nothing"""

    # These statements check to make sure that all arguments were entered correctly.
    assert isinstance(csv_file, (str, unicode))
    assert isinstance(filename, (str, unicode))
    assert isinstance(D, int)
    assert isinstance(T, int)
    assert isinstance(L, int)
    assert isinstance(debug, bool)
    assert isinstance(zero_loss, bool)

    # Record the time that DP starts
    start_time = time.clock()

    # Get everything we need from DP
    edge_species_tree, edge_gene_tree, graph, mpr_count = DP.reconcile(filename, D, T, L)

    # And record the amount of time DP took
    DP_time_taken = time.clock() - start_time

    print "Reconciliation Complete in \033[33m\033[1m{0} seconds\033[0m".format(DP_time_taken)

    # Record the time that this code starts
    start_time = time.clock()

    # The gene tree needs to be in node format, not edge format, so we find that now.
    # (This also puts the gene_tree into postorder, as an ordered dict)
    gene_tree, gene_tree_root, gene_node_count = reformat_tree(edge_gene_tree, "pTop")

    # The DTL reconciliation graph as provided by DP has some extraneous numbers and may have loss events on the root
    # node. We remove those here.
    clean_graph(graph, gene_tree_root)

    # Now we draw the rest of the owl
    diameter = diameter_algorithm(edge_species_tree, gene_tree, gene_tree_root, graph, debug, zero_loss)

    # And record how long it took to compute the diameter.
    diameter_time_taken = time.clock() - start_time

    print "The diameter of the given reconciliation graph is \033[33m\033[1m{0}\033[0m".format(diameter)

    # Timing data is inaccurate in debug mode (print statements take too long), so we only give it to the user in non-
    # debug mode.
    if not debug:
        print "Done in \033[33m\033[1m{0} seconds\033[0m".format(diameter_time_taken)

    # Now, we write our results to a csv file.
    file_exists = os.path.isfile(csv_file)

    with open(csv_file, 'a') as output_file:
        writer = csv.writer(output_file)

        # Write the headers if we need to.
        if not file_exists:
            writer.writerow(["File Name", "MPR Count", "Diameter Found", "Gene Node Count",
                             "DP Computation Time", "Diameter Computation Time", "Date"])

        writer.writerow([filename, mpr_count, diameter, gene_node_count,
                         DP_time_taken, diameter_time_taken, time.strftime("%c")])

    # And we're done.
    return

def c(file_name="example", D=0, T=0, L=0, test_log="TestLog.csv", Debug=False, zero_loss=False):
    """This method is identical to calculate_diameter_from_file, but with a shorter name that's easier to type"""
    calculate_diameter_from_file(file_name, D, T, L, test_log, Debug, zero_loss)


def t(file_name="example"):
    """A function to call calculate_diameter with some testing values, because typing that name fully is slower overall
     than writing this function (and this docstring)"""
    calculate_diameter_from_file(file_name, 0, 0, 0, "TestLog.csv", True)
    print "Expected value: 9"


def t2(file_name="example"):
    calculate_diameter_from_file(file_name, 1, 4, 1, "TestLog.csv", True)
    print "Expected value: 7"