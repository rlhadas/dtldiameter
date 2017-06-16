import DP

# TODO: Find out what Jordan's original algorithm for finding the intersections between two paths was, and see if mine is better
# TODO: Find out what Group(u) means computationally, and whether that is better than what I have done

def ReformatTree(Tree,Root):
    """A recursive function that changes the format of a (species or gene) tree from edge to vertex, for example:
    ('A','B'): ('A','B',('B',C1),('B',C2)) would become 'B':(C1,C2). It returns the tree and the root."""

    NewRoot = Root[1] if isinstance(Root, tuple) else Tree[Root][1] # This line catches the "xTop" handle

    Child1 = Tree[Root][2][1] if Tree[Root][2] is not None else None # These lines handle the leaves, where
    Child2 = Tree[Root][3][1] if Tree[Root][3] is not None else None # there is None in the place of a tuple

    NewTree = {NewRoot:(Child1, Child2)}
    if Child1 is not None: # If this node has children, then we need to add their children's subtrees to the dict
        Child1Tree, _ = ReformatTree(Tree, Tree[Root][2])
        NewTree.update(Child1Tree)
    if Child2 is not None:
        Child2Tree, _ = ReformatTree(Tree, Tree[Root][3])
        NewTree.update(Child2Tree)

    return NewTree, NewRoot


def FindValidPaths(Root, PreviousValues, Tree):
    """A recursive algorithm to find all of the valid paths through a binary tree, as well as a dict containing the
    nodes in each path."""

    # The tree is represented as a dict, with values having the structure (previousNode,thisNode,childNode1,childNode2)
    # With the child nodes values being the keys to those "NodeSet"s, as I've called them here

    NextValues = PreviousValues + [Root]  # This list contains every node we visited to get here from the root
    Paths = []
    PathNodes = {}
    for i in range(0, len(NextValues)):
        SourceNode = NextValues[i]
        NewPath = (SourceNode, Root)  # This becomes every path that ends in this value, including A->A
        Paths += [NewPath]
        PathNodes[NewPath] = NextValues[i:]  # We need to add every node after and including the SourceNode

    Child1 = Tree[Root][0]
    Child2 = Tree[Root][1]

    if not Child1 == None:  # Then this Node is not a leaf Node, so we need to add this Node's children
        Child1Paths, Child1PathNodes = FindValidPaths(Child1, NextValues, Tree)
        Child2Paths, Child2PathNodes = FindValidPaths(Child2, NextValues, Tree)
        Paths += Child1Paths + Child2Paths
        Child1PathNodes.update(Child2PathNodes)
        PathNodes.update(Child1PathNodes)
                                # Otherwise, we have reached the end of the tree (the base case)
    return Paths, PathNodes


def ComputePathSymmetricSetDifferenceTable(SpeciesTree, SpeciesTreeRoot):
    """Computes the table containing the number of nodes in the symmetric set difference between any two paths on the
     species tree sTree. This is used in assigning a score to two paths' losses."""

    PathList, PathNodes = FindValidPaths(SpeciesTreeRoot, [], SpeciesTree)

    # Note: If you wish to modify the code to assign different scores for each gene node, modifying this function
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
        SSD[PathA] = {}
        for PathB in NonTrivialPathList:
            BNodes = frozenset(PathNodes[PathB])
            SSD[PathA][PathB] = len(ANodes.symmetric_difference(BNodes))

    return SSD, PathList, NonTrivialPathList


def BuildExitDicts(Graph):
    """Builds a dict containing lists of exit event nodes, each list keyed by the gene node."""
    # The graph returned by DP is keyed by mapping nodes, so to get the gene node we take the first element of the key.
    ExitEventDict = {} # Every gene node's list of exit events
    ExitMappingDict = {} # Every gene node's list of mapping nodes that have exit events
    for MappingNode in Graph:
        for Event in Graph[MappingNode]:
            if isinstance(Event,tuple) and Event[0] != 'L': # Ignore Loss Events and Scores
                GeneNode = MappingNode[0]
                if GeneNode not in ExitEventDict:
                    ExitEventDict[GeneNode] = []
                ExitEventDict[GeneNode].append(Event)

                if GeneNode not in ExitMappingDict:
                    ExitMappingDict[GeneNode] = []
                ExitMappingDict[GeneNode] += [MappingNode]

    return ExitEventDict, ExitMappingDict


def ComputeTrivialExitEventTable(u, ExitEvent):
    """This function computes and stores the score of the exit event on a leaf node 'u' of the gene tree.
    As this event will always be a C event that is shared by all nodes, this value will always be 0."""
    ExitEvent[u] = {}
    ExitEvent[u][('C', (None, None), (None, None))] = {}
    ExitEvent[u][('C', (None, None), (None, None))][('C', (None, None), (None, None))] = 0


def ComputeExitEventTable(u, ExitEvent, EnterMapping, ExitEventList):
    """This function computes and stores the score of the exit event on a non-leaf node 'u' of the gene tree."""
    #TODO: Test this function
    for E1 in ExitEventList[u]:
        Child1 = E1[1][0]
        Child2 = E1[2][1]
        uB = E1[1]
        uC = E1[2]
        ExitEvent[u][E1] = {}
        for E2 in ExitEventList[u]:
            #We need to account for the case that the children of u are in opposite order between the two events
            if Child1 == E2[1][0]:
                uE = E2[1]
                uF = E2[2]
            else:
                uE = E2[2]
                uF = E2[1]
            print "E1:" + str(E1)
            print "E2:" + str(E2)
            ExitEvent[u][E1][E2] = EnterMapping[Child1][uB][uE] \
                                   + EnterMapping[Child2][uC][uF] \
                                   + 1 if E1 != E2 else 0



def ComputeExitMappingTable(u, ExitMapping, ExitMappingNodeDict):
    """This function computes and stores the maximum possible score of the exit from gene node u"""

    # TODO: Replace filler in function

    uMappingNodes = ExitMappingNodeDict[u]
    ExitMapping[u] = {}
    for uA in uMappingNodes:
        ExitMapping[u][uA] = {}
        for uB in uMappingNodes:
            ExitMapping[u][uA][uB] = 2 # TODO: This is filler, replace with actual algorithm


def ComputeEnterMappingTable(u, EnterMapping, MappingNodeList):
    """This function computes the maximum possible score of each pair of mapping nodes for gene node u, and stores each
    one into the EnterMapping table for u."""

    # TODO: Replace filler in function

    uMappingNodes = [] # Make a new list that has only the mapping nodes that contain u
    for node in MappingNodeList:
        if node[0] == u:
            uMappingNodes.append(node)
    # Todo: create LossReachable list

    EnterMapping[u] = {}
    for uA in uMappingNodes:
        EnterMapping[u][uA] = {}
        for uB in uMappingNodes:
            EnterMapping[u][uA][uB] = 2 # TODO: This is filler, replace it with the right algorithm.

def PrintTableNicely(Table, Deliminator, Name="\t"):
    """Takes a table (a 2D dict keyed with tuples) and prints a nicely formatted table. Used for debugging."""
    line = "{0}\t".format(Name)
    for Column in Table:
        line += "({0}{1}{2})\t".format(str(Column[0]), Deliminator, str(Column[1]))
    print line
    for Row in Table:
        line = "\t({0}{1}{2}):\t".format(str(Row[0]), Deliminator, str(Row[1]))
        for Column in Table:
            line += str(Table[Column][Row]) + "\t"
        print line

def CalculateDiameter(filename, D, T, L):
    """This function computes the diameter of space of MPRs in a DTL reconciliation problem,
    as measured by the symmetric set distance between the two."""

    # TODO: Add a way to run the CalculateDiameter function on a previously-found reconciliation.

    EdgeSpeciesTree, EdgeGeneTree, Graph = DP.reconcile(filename, D, T, L)
    print "Reconciliation Complete"

    SpeciesTree, SpeciesTreeRoot = ReformatTree(EdgeSpeciesTree, "hTop")
    GeneTree, GeneTreeRoot = ReformatTree(EdgeGeneTree, "pTop")

    for key in Graph:  # The events in the event graph are stored as lists which cannot be used as dict keys.
        for i in range(0, len(Graph[key])):
            if isinstance(Graph[key][i], list):  # Ignore the "scores" for now
                Graph[key][i] = tuple(Graph[key][i])

    PathSymmetricSetDifference = {} # A dict containing the SSD count for each pair of species node paths.

    PathList = []   # A list of all of the valid paths between two species nodes (format is (src,dest))

    NonTrivialPathList = []  # A subset of the elements in PathList where src != dest

    # The key format for the next three dicts are as follows: (u, x, y), where u is a node on the gene tree, and
    # x and y are either event nodes (represented as E1 and E2) or mapping nodes (represented as uA and uB).
    # TODO: Clear entries of these dicts when they are no longer needed

    ExitEvent = {}  # A dict containing the largest number of event nodes that each pair of reconciliation subtrees
                    # rooted at E1 and E2 can differ for, whereE1 and E2 are exit-event nodes in Group(u).

    ExitMapping = {}    # A dict containing the largest number of event nodes that each pair of reconciliation subtrees
                        # rooted at uA and uB can have in common, where uA and uB go immediately to an exit-event.

    EnterMapping = {}   # A dict containing the largest number of event nodes that each pair of reconciliation subtrees
                        # rooted at uA and uB can have in common.

    ExitEventDict = {}  # This dict contains each exit (non-loss) event in the reconciliation graph, keyed by gene node.

    MappingNodeList = Graph.keys() # This list contains each mapping node in the reconciliation graph


    PathSymmetricSetDifference, PathList, NonTrivialPathList = \
        ComputePathSymmetricSetDifferenceTable(SpeciesTree, SpeciesTreeRoot)


    ExitEventDict, ExitMappingNodeDict = BuildExitDicts(Graph)
    PrintTableNicely(PathSymmetricSetDifference,"->","[[SSD]]:")
    print ExitMappingNodeDict
    PostorderGeneVertices = GeneTree.keys() # TODO: get this sorted in postorder
    PostorderGeneVertices.reverse() # TODO: check to see if this is always in postorder (probably not)
    for u in PostorderGeneVertices:
        EnterMapping[u] = {}
        ExitMapping[u] = {}
        ExitEvent[u] = {}

    for u in PostorderGeneVertices:
        if GeneTree[u][0] is not None:  # Then u is not a leaf node
            assert EnterMapping[GeneTree[u][0]] != {} and \
                EnterMapping[GeneTree[u][0]] != {}
            ComputeExitEventTable(u, ExitEvent, EnterMapping, ExitEventDict)
        else:  # u IS a leaf node, and therefore its ExitEvent table ExitEvent[u] is trivial
            ComputeTrivialExitEventTable(u, ExitEvent)
        ComputeExitMappingTable(u, ExitMapping, ExitMappingNodeDict)
        ComputeEnterMappingTable(u, EnterMapping, MappingNodeList)

        PrintTableNicely(ExitEvent[u], ", ", "ExitEvent({0})\t\t".format(u))
        PrintTableNicely(ExitMapping[u], ", ", "ExitMapping({0})".format(u))
        PrintTableNicely(EnterMapping[u], ", ", "EnterMapping({0})".format(u))




def tst():
    CalculateDiameter("example", 0, 0, 0)