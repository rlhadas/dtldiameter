import DP

def FindValidPaths(Key, PreviousValues, Tree):
    """A recursive algorithm to find all of the valid paths through a binary tree, as well as a dict containing the
    nodes in each path."""

    #The tree is represented as a dict, with values having the structure (previousNode,thisNode,childNode1,childNode2)
    #With the child nodes values being the keys to those "NodeSet"s, as I've called them here

    NodeSet = Tree[Key]
    Node = NodeSet[1]
    NextValues = PreviousValues + [Node]  # This list contains every node we visited to get here from the root
    Paths = []
    PathNodes = {}
    for i in range(0, len(NextValues)):
        SourceNode = NextValues[i]
        NewPath = (SourceNode, Node)  # This becomes every path that ends in this value, including A->A
        Paths += [NewPath]
        PathNodes[NewPath] = NextValues[i:]  # We need to add every node after and including the SourceNode

    if not NodeSet[2] == None:  # Then this Node is not a leaf Node, so we need to add this Node's children
        Child1Paths, Child1PathNodes = FindValidPaths(NodeSet[2], NextValues, Tree)
        Child2Paths, Child2PathNodes = FindValidPaths(NodeSet[3], NextValues, Tree)
        Paths += Child1Paths + Child2Paths
        Child1PathNodes.update(Child2PathNodes)
        PathNodes.update(Child1PathNodes)
                                # Otherwise, we have reached the end of the tree (the base case)
    return Paths, PathNodes



def ComputePathSymmetricSetDifferenceTable(sTree):
    """Computes the table containing the number of nodes in the symmetric set difference between any two paths on the
     species tree sTree. This is used in assigning a score to two paths' losses."""

    PathList, PathNodes = FindValidPaths('hTop', [], sTree) # The handle of sTree is 'hTop', as given by newickFormatReader.py

    # Note: If you wish to modify the code to assign different costs for each gene node, modifying this function
    # to provide the actual lists of nodes in the SSD might be a good place to start

    SSD = {} #This is the Symmetric Set Difference table we will be returning

    # The algorithm described in Jordan's writeup to find the intersection did not make sense to me, and I don't know
    # whether it applies to the SSD

    NonTrivialPathList = list(PathList) # Strip the trivial paths because they will not have loss events
    for Path in NonTrivialPathList:
        if Path[0] == Path[1]:
            NonTrivialPathList.remove(Path)

    for PathA in NonTrivialPathList:
        ANodes = frozenset(PathNodes[PathA])
        for PathB in NonTrivialPathList:
            BNodes = frozenset(PathNodes[PathB])
            SSD[(PathA, PathB)] = len(ANodes.symmetric_difference(BNodes))

    return SSD, PathList, NonTrivialPathList


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


def PrintPathTableNicely(PathList,PathTable):
    line = "\t"
    for Px in PathList:
        line += str(Px[0]) + "->" + str(Px[1]) + "\t"
    print line
    for Py in PathList:
        line = str(Py[0]) + "->" + str(Py[1]) + ":\t"
        for Px in PathList:
            line += str(PathTable[(Px,Py)]) +"\t"
        print line

def CalculateDiameter(filename, D, T, L):
    """This function computes the diameter of space of MPRs in a DTL reconciliation problem,
    as measured by the symmetric set distance between the two."""

    # TODO: Add a way to run the CalculateDiameter function on a previously-found reconciliation.

    SpeciesTree, RecTree = DP.reconcile(filename, D, T, L)
    print "Reconciliation Complete"
    PathSymmetricSetDifference = {} # A dict containing the SSD count for each pair of species node paths.

    PathList = []   # A list of all of the valid paths between two species nodes (format is (src,dest))

    NonTrivialPathList = []  # A subset of the elements in PathList where src != dest

    # The key format for the next three dicts are as follows: (u, x, y), where u is a node on the gene tree, and
    # x and y are either event nodes (represented as E1 and E2) or mapping nodes (represented as uA and uB).

    ExitEvent = {}  # A dict containing the largest number of event nodes that each pair of reconciliation subtrees
                    # rooted at E1 and E2 can differ for, whereE1 and E2 are exit-event nodes in Group(u).

    ExitMapping = {}    # A dict containing the largest number of event nodes that each pair of reconciliation subtrees
                        # rooted at uA and uB can have in common, where uA and uB go immediately to an exit-event.

    EnterMapping = {}   # A dict containing the largest number of event nodes that each pair of reconciliation subtrees
                        # rooted at uA and uB can have in common.

    PathSymmetricSetDifference, PathList, NonTrivialPathList = ComputePathSymmetricSetDifferenceTable(SpeciesTree)
    print PathList

