import DP


def ComputePathSymmetricSetDifferenceTable():
    """Computes the table containing the number of nodes in the symmetric set difference between any two paths on the
     species tree."""
    # TODO: everything
    # Note: If you wish to modify the code to assign different costs for each gene node, modifying this function
    # to provide the actual lists of nodes in the SSD might be a good place to start
    return {}, []


def ComputeTrivialExitEventTable(u, ExitEvent):
    """This function computes and stores the cost of the exit event on a leaf node 'u' of the gene tree.
    As this event will always be a C event that is shared by all nodes, this value will always be 0."""

    #TODO: Add function body



def ComputeExitEventTable(u):
    """This function computes and stores the cost of the exit event on a non-leaf node 'u' of the gene tree."""

    # TODO: Add function body


def ComputeExitMappingTable(u):
    """This function computes and stores the maximum possible cost of the exit from gene node u"""

    # TODO: Add function body


def ComputeEnterMappingTable(u):
    """This function computes the maximum possible cost of each pair of mapping nodes for gene node u, and stores each
    one into the EnterMapping table for u."""

    # TODO: Add function body


def CalculateDiameter(filename, D, T, L):
    """This function computes the diameter of space of MPRs in a DTL reconciliation problem,
    as measured by the symmetric set distance between the two."""

    Tree = DP.reconcile(filename, D, T, L)

    PathSymmetricSetDifference = {} # A dict containing the SSD for each non-trivial pair of species nodes.

    PathList = []   # A list of all of the valid paths between two species nodes

    # The key format for the next three dicts are as follows: (u, x, y), where u is a node on the gene tree, and
    # x and y are either event nodes (represented as E1 and E2) or mapping nodes (represented as uA and uB).

    ExitEvent = {}  # A dict containing the largest number of event nodes that each pair of reconciliation subtrees
                    # rooted at E1 and E2 can differ for, whereE1 and E2 are exit-event nodes in Group(u).

    ExitMapping = {}    # A dict containing the largest number of event nodes that each pair of reconciliation subtrees
                        # rooted at uA and uB can have in common, where uA and uB go immediately to an exit-event.

    EnterMapping = {}   # A dict containing the largest number of event nodes that each pair of reconciliation subtrees
                        # rooted at uA and uB can have in common.

    PathSymmetricSetDifference, PathList = ComputePathSymmetricSetDifferenceTable()

