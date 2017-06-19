import DP


# TODO: Find out what Jordan's original algorithm for finding the intersections between two paths was, and see if mine is better
# TODO: Find out what Group(u) means computationally, and whether that is better than what I have done

def reformat_tree(tree, root):
    """A recursive function that changes the format of a (species or gene) tree from edge to vertex, for example:
    ('A','B'): ('A','B',('B',C1),('B',C2)) would become 'B':(C1,C2). It returns the tree and the root."""

    new_root = root[1] if isinstance(root, tuple) else tree[root][1]  # This line catches the "xTop" handle

    child1 = tree[root][2][1] if tree[root][2] is not None else None  # These lines handle the leaves, where
    child2 = tree[root][3][1] if tree[root][3] is not None else None  # there is None in the place of a tuple

    new_tree = {new_root: (child1, child2)}
    if child1 is not None:  # If this node has children, then we need to add their children's subtrees to the dict
        Child1Tree, _ = reformat_tree(tree, tree[root][2])
        new_tree.update(Child1Tree)
    if child2 is not None:
        Child2Tree, _ = reformat_tree(tree, tree[root][3])
        new_tree.update(Child2Tree)

    return new_tree, new_root


def find_valid_paths(root, previous_values, tree):
    """A recursive algorithm to find all of the valid paths through a binary tree, as well as a dict containing the
    nodes in each path."""

    next_values = previous_values + [root]  # This list contains every node we visited to get here from the root
    paths = []
    path_nodes = {}
    for i in range(0, len(next_values)):
        source_node = next_values[i]
        new_path = (source_node, root)  # This becomes every path that ends in this value, including A->A
        paths += [new_path]
        path_nodes[new_path] = next_values[i:]  # We need to add every node after and including the source_node

    child1 = tree[root][0]
    child2 = tree[root][1]

    if child1 is not None:  # Then this Node is not a leaf Node, so we need to add this Node's children
        child1_paths, child1_path_nodes = find_valid_paths(child1, next_values, tree)
        child2_paths, child2_path_nodes = find_valid_paths(child2, next_values, tree)
        paths += child1_paths + child2_paths
        child1_path_nodes.update(child2_path_nodes)
        path_nodes.update(child1_path_nodes)
        # Otherwise, we have reached the end of the tree (the base case)
    return paths, path_nodes


def compute_path_symmetric_set_difference_table(species_tree, species_tree_root):
    """Computes the table containing the number of nodes in the symmetric set difference between any two paths on the
     species tree sTree. This is used in assigning a score to two paths' losses."""

    path_list, path_nodes = find_valid_paths(species_tree_root, [], species_tree)

    # Note: If you wish to modify the code to assign different scores for each gene node, modifying this function
    # to provide the actual lists of nodes in the SSD might be a good place to start

    ssd = {}  # This is the Symmetric Set Difference table we will be returning

    # The algorithm described in Jordan's writeup to find the intersection did not make sense to me, and I don't know
    # whether it applies to the SSD

    non_trivial_path_list = list(path_list)  # Strip the trivial paths because they will not have loss events
    for path in non_trivial_path_list:
        if path[0] == path[1]:
            non_trivial_path_list.remove(path)

    for path_a in non_trivial_path_list:
        a_nodes = frozenset(path_nodes[path_a])
        ssd[path_a] = {}
        for path_b in non_trivial_path_list:
            b_nodes = frozenset(path_nodes[path_b])
            ssd[path_a][path_b] = len(a_nodes.symmetric_difference(b_nodes))

    return ssd, path_list, non_trivial_path_list


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


def compute_trivial_exit_event_table(u, exit_event):
    """This function computes and stores the score of the exit event on a leaf node 'u' of the gene tree.
    As this event will always be a C event that is shared by all nodes, this value will always be 0."""
    exit_event[u] = {}
    exit_event[u][('C', (None, None), (None, None))] = {}
    exit_event[u][('C', (None, None), (None, None))][('C', (None, None), (None, None))] = 0


def compute_exit_event_table(u, exit_event, enter_mapping, exit_event_list):
    """This function computes and stores the score of the exit event on a non-leaf node 'u' of the gene tree."""

    # TODO: Make sure this function provides accurate values

    for e1 in exit_event_list[u]:
        child1 = e1[1][0]
        child2 = e1[2][0]
        # B and C are the species nodes of the two mapping nodes of e1
        B = e1[1][1]
        C = e1[2][1]
        exit_event[u][e1] = {}
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
            exit_event[u][e1][e2] = enter_mapping[child1][uB][uE] \
                                    + enter_mapping[child2][uC][uF] \
                                    + (1 if e1 != e2 else 0)


def compute_exit_mapping_table(u, exit_mapping, exit_event, exit_mapping_node_dict, exit_event_graph):
    """This function computes and stores the maximum possible score of the exit from gene node u"""

    # TODO: Make sure this function provides accurate values

    u_mapping_nodes = exit_mapping_node_dict[u]
    exit_mapping[u] = {}
    for uA in u_mapping_nodes:
        exit_mapping[u][uA] = {}
        for uB in u_mapping_nodes:
            max_value = 0
            for E1 in exit_event_graph[uA]:
                for E2 in exit_event_graph[uB]:
                    max_value = max(exit_event[u][E1][E2], max_value)
            exit_mapping[u][uA][uB] = max_value


def create_loss_reachable(graph, root_mapping_node):
    """A recursive function to create the list of the mapping nodes that can be reached through loss events.
    (The base case is when you can't get to a loss node from a mapping node)"""
    loss_reachable = []

    loss_events = filter(lambda e: e[0] == 'L', graph[root_mapping_node])  # Grab all loss events
    for Event in loss_events:
        print Event
        loss_reachable += create_loss_reachable(graph, Event[1])
    return loss_reachable + [root_mapping_node]  # TODO: Add this node (and make sure the overall root node is not added)


def compute_enter_mapping_table(u, enter_mapping, mapping_node_list, graph, ssd):
    """This function computes the maximum possible score of each pair of mapping nodes for gene node u, and stores each
    one into the enter_mapping table for u."""

    # TODO: Replace filler in function

    u_mapping_nodes = []  # Make a new list that has only the mapping nodes that contain u
    for node in mapping_node_list:
        if node[0] == u:
            u_mapping_nodes.append(node)

    loss_reachable = {}  # Every mapping node (in u)'s list of mapping nodes that can be reached through loss events

    for mapping_node in u_mapping_nodes:
        loss_reachable[mapping_node] = create_loss_reachable(graph, mapping_node)
        loss_reachable[mapping_node].remove(mapping_node)  # You can't reach the root mapping node from itself
        print loss_reachable[mapping_node]

    enter_mapping[u] = {}
    for uA in u_mapping_nodes:
        enter_mapping[u][uA] = {}
        for uB in u_mapping_nodes:
            enter_mapping[u][uA][uB] = 2  # TODO: This is filler, replace it with the right algorithm.


def event_to_string(event):
    return "{0}:{1}{2} {3}{4}".format(str(event[0]), str(event[1][0]), str(event[1][1]),
                                      str(event[2][0]), str(event[2][1]))


def print_table_nicely(table, deliminator, name="\t", is_event=False):
    """Takes a table (a 2D dict keyed with tuples) and prints a nicely formatted table. Used for debugging."""
    print ""
    line = "{0}".format(name)
    for column in table:
        if is_event:
            line += "\t{0}".format(event_to_string(column))
        else:
            line += "\t{0}{1}{2}".format(str(column[0]), deliminator, str(column[1]))
    print line
    for row in table:
        if is_event:
            line = "\t{0}\t".format(event_to_string(row))
        else:
            line = "\t{0}{1}{2}\t".format(str(row[0]), deliminator, str(row[1]))
        for column in table:
            line += str(table[column][row]) + "\t"
        print line


def sanitize_graph(graph):
    """Cleans up the graph created by DP.py by turning events into tuples and removing scores from events.
    This allows us to use the events as dictionary keys."""

    for key in graph:
        # Get rid of all of the random numbers in the graph
        graph[key] = filter(lambda e: not isinstance(e, int) and not isinstance(e, float), graph[key])
        # The events in the event graph are stored as lists which cannot be used as dict keys. Let's fix that.
        for i in range(0, len(graph[key])):
            if isinstance(graph[key][i], list):
                event = graph[key][i]
                event = filter(lambda x: not isinstance(x, float), event)  # Remove more numbers
                graph[key][i] = tuple(event)


def calculate_diameter(filename, D, T, L):
    """This function computes the diameter of space of MPRs in a DTL reconciliation problem,
    as measured by the symmetric set distance between the events of the two reconciliations of the pair
     that has the highest such difference."""

    # TODO: Add a way to run the calculate_diameter function on a previously-found reconciliation.

    edge_species_tree, edge_gene_tree, graph = DP.reconcile(filename, D, T, L)
    print "Reconciliation Complete"

    species_tree, species_tree_root = reformat_tree(edge_species_tree, "hTop")
    gene_tree, gene_tree_root = reformat_tree(edge_gene_tree, "pTop")

    sanitize_graph(graph)

    path_symmetric_set_difference = {}  # A dict containing the SSD count for each pair of species node paths.

    path_list = []  # A list of all of the valid paths between two species nodes (format is (src,dest))

    non_trivial_path_list = []  # A subset of the elements in path_list where src != dest

    # The key format for the next three dicts are as follows: (u, x, y), where u is a node on the gene tree, and
    # x and y are either event nodes (represented as E1 and E2) or mapping nodes (represented as uA and uB).
    # TODO: Clear entries of these dicts when they are no longer needed

    exit_event_scores = {}  # A dict containing the largest number of event nodes that each pair of reconciliation subtrees
    # rooted at E1 and E2 can differ for, whereE1 and E2 are exit-event nodes in Group(u).

    exit_mapping_scores = {}  # A dict containing the largest number of event nodes that each pair of reconciliation subtrees
    # rooted at uA and uB can have in common, where uA and uB go immediately to an exit-event.

    enter_mapping_scores = {}  # A dict containing the largest number of event nodes that each pair of reconciliation subtrees
    # rooted at uA and uB can have in common.

    exit_event_dict = {}  # This dict contains each exit (non-loss) event in the reconciliation graph, keyed by gene node.

    exit_event_graph = {}  # This dict contains each exit (non-loss) event in the reconciliation graph, keyed by mapping node.

    mapping_node_list = graph.keys()  # This list contains each mapping node in the reconciliation graph

    path_symmetric_set_difference, path_list, non_trivial_path_list = \
        compute_path_symmetric_set_difference_table(species_tree, species_tree_root)

    exit_event_graph, exit_event_dict, exit_mapping_node_dict = build_exit_dicts(graph)
    print_table_nicely(path_symmetric_set_difference, "->", "[[SSD]]:")

    postorder_gene_vertices = gene_tree.keys()
    postorder_gene_vertices.reverse()  # TODO: check to see if this is always in postorder (probably not)
    for u in postorder_gene_vertices:
        enter_mapping_scores[u] = {}
        exit_mapping_scores[u] = {}
        exit_event_scores[u] = {}

    for u in postorder_gene_vertices:
        if gene_tree[u][0] is not None:  # Then u is not a leaf node
            assert enter_mapping_scores[gene_tree[u][0]] != {} and \
                   enter_mapping_scores[gene_tree[u][0]] != {}
            compute_exit_event_table(u, exit_event_scores, enter_mapping_scores, exit_event_dict)
        else:  # u IS a leaf node, and therefore its exit_event_scores table exit_event_scores[u] is trivial
            compute_trivial_exit_event_table(u, exit_event_scores)
        compute_exit_mapping_table(u, exit_mapping_scores, exit_event_scores, exit_mapping_node_dict, exit_event_graph)
        compute_enter_mapping_table(u, enter_mapping_scores, mapping_node_list, graph, path_symmetric_set_difference)

        print_table_nicely(exit_event_scores[u], ", ", "exit_event_scores({0})".format(u), True)
        print_table_nicely(exit_mapping_scores[u], "", "exit_mapping_scores({0})".format(u))
        print_table_nicely(enter_mapping_scores[u], "", "enter_mapping_scores({0})".format(u))


def t():
    """A function to call calculate_diameter with some testing values, because typing that name fully is slower overall
     than writing this function (and this docstring)"""
    return calculate_diameter("example", 0, 0, 0)
