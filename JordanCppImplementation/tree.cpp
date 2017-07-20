// Author: Jordan Haack

#include <string>
#include <list>
#include <utility>
#include <unordered_map>
#include <vector>
#include <queue>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdexcept>

// Easy int to string conversions
#include <sstream>
#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

using namespace std;

//-------|---------|---------|---------|---------|---------|---------|---------|
/// A vertex for a tree
class Vertex {
public:

  	Vertex *left, *right;
	string label;
    Vertex* parent = nullptr;
    
    size_t index; // indices are assigned in post order
    int height;

    Vertex(string qlabel) 
        : left(nullptr), right(nullptr), label(qlabel) { }

    Vertex(Vertex* qleft, Vertex* qright, string qlabel) 
        : left(qleft), right(qright), label(qlabel) { }

    ~Vertex() {
        if (left) delete left;
        if (right) delete right;
    }

    bool isLeaf() {
        return left == nullptr;
    }

    bool isRoot() {
        return parent == nullptr;
    }

    void print(string& s) {
        if (left) {
            s += "(";
            left->print(s);
            s += ",";
            right->print(s);
            s += ")" + label;
        } else {
            s += label;
        }
    }

    /// A tree helper function, to learn info about the tree,
    /// as well as assign indexes and finding labels, height, size, etc.
    void populateInfo(size_t& size, vector<Vertex*>& v, 
                        unordered_map<string, Vertex*>& m) {
        if (left) {
            left->populateInfo(size, v, m);
            right->populateInfo(size, v, m);
            left->parent = this;
            right->parent = this;
            height = 1 + max<int>(left->height, right->height);
        } else {
            height = 0;
        }
        index = size;
        size++;
        v.push_back(this);
        m[label] = this;
    }

    void preOrder(vector<Vertex*>& L) {
        L.push_back(this);
        if (left) {
            left->preOrder(L);
            right->preOrder(L);
        }
    }

    void postOrder(vector<Vertex*>& L) {
        if (left) {
            left->postOrder(L);
            right->postOrder(L);
        }
        L.push_back(this);
    }

    size_t inOrderIndex; // Used for drawing a tree
    void inOrder(vector<Vertex*>& L) {
        if (left) 
            left->inOrder(L);
        L.push_back(this);
        if (right)
            right->inOrder(L);
    }

    // Used as storage for drawing trees.
    int posX, posY;
};


//-------|---------|---------|---------|---------|---------|---------|---------|

enum class Rel {
    EQUAL = 2,
    ANCESTOR = 1,
    DESCENDENT = -1,
    INCOMPARABLE = 0,
};

/// A tree class, made up of vertices. Used to represent the
/// host and parasite trees
class Tree {
public:
	Vertex* root;

    size_t size;
    int height;
    vector<Vertex*> vertexList; // Post order list
    vector<Vertex*> preOrderVertexList; // Pre order list
    vector<Vertex*> inOrderVertexList; // in order list
    unordered_map<string, Vertex*> labelMap; // Used to compute the tipMapping

    Tree(Vertex* qroot) : root(qroot) { 
        populateInfo();
        computeAncestorTable();
    }

    ~Tree() {
        if (root) delete root;
    }

    string printToString() {
        string s;
        root->print(s);
        return s;
    }

    /// Updates information about the tree such as size, labels,
    /// and pre/post order traversals.
    void populateInfo() {
        size = 0;
        vertexList = vector<Vertex*>();
        labelMap = unordered_map<string, Vertex*>();
        root->populateInfo(size, vertexList, labelMap);
        height = root->height;
        preOrderVertexList = preOrderTraversal();
        inOrderVertexList = inOrderTraversal();
        int ind = 0;
        for (Vertex* v: inOrderVertexList)
            v->inOrderIndex = ind++;
    }

    /// PreOrder. O(n) time
    vector<Vertex*> preOrderTraversal() {
        vector<Vertex*> L;
        root->preOrder(L);
        return L;
    }

    /// PostOrder. O(n) time
    vector<Vertex*> postOrderTraversal() {
        vector<Vertex*> L;
        root->postOrder(L);
        return L;
    }

    /// InOrder. O(n) time
    vector<Vertex*> inOrderTraversal() {
        vector<Vertex*> L;
        root->inOrder(L);
        return L;
    }

    /// The ancestor table takes the index of two nodes and 
    /// allows you to see their relationship on the tree.
    /// It's kinda hacky but it works
    /// O(n^2) time
    vector<vector<Rel>> ancestorTable;
    void computeAncestorTable() {
        ancestorTable = vector<vector<Rel>>(size, vector<Rel>(size));
        for (size_t a = size-1; a < size; --a) { // pre order, note a is unsigned
            for (size_t d = size-1; d < size; --d) { // pre order
                if (vertexList[a] == vertexList[d]) 
                    ancestorTable[a][d] = Rel::EQUAL;
                else if (ancestorTable[d][a] == Rel::ANCESTOR)
                    ancestorTable[a][d] = Rel::DESCENDENT;
                else if (vertexList[a] == vertexList[d]->parent)
                    ancestorTable[a][d] = Rel::ANCESTOR;
                else if (ancestorTable[a][vertexList[d]->parent->index] == Rel::ANCESTOR)
                    ancestorTable[a][d] = Rel::ANCESTOR;
                else
                    ancestorTable[a][d] = Rel::INCOMPARABLE;
            }
        }
    }

};

//-------|---------|---------|---------|---------|---------|---------|---------|
/// All of the possible event node types are here
enum class EventType {
    C, D, T, L, S
};

//-------|---------|---------|---------|---------|---------|---------|---------|
/// The nodes found on a reconciliation graph
class Node{
public:
    class MappingNode;
    class EventNode;

  /// Event node class for storing event nodes on a reconciliaiton graph
  class EventNode {
  public:
    EventType type;
    MappingNode *left, *right, *parent;

    EventNode(EventType t, MappingNode* l, MappingNode* r, MappingNode* p) 
    : type(t), left(l), right(r), parent(p) { }

    // Utilities
    int dpScore;
    double bottomUpScore; // Used in frequency calculation
    double numberOfTrees; // Number or trees appeared in
    float frequency;

    string str() {
        string s = "";
        switch (type) {
            case EventType::C: s += "Cont"; break;
            case EventType::D: s += "Dupl"; break;
            case EventType::T: s += "Trns"; break;
            case EventType::L: s += "Loss"; break;
            case EventType::S: s += "Spec"; break;
        }

        s += " from " + parent->str();
        if (left != nullptr) s += " to "+ left->str();
        if (right != nullptr) s += " and " + right->str();
        return s;
    }

    string cstr() {  // Condensed string
        string s = "";
        switch (type) {
            case EventType::C: s += "C"; break;
            case EventType::D: s += "D"; break;
            case EventType::T: s += "T"; break;
            case EventType::L: s += "L"; break;
            case EventType::S: s += "S"; break;
        }
        if (left != nullptr) s += ":"+ left->cstr();
        if (right != nullptr) s += ":" + right->cstr();
        return s;
    }

    // Used as storage for drawing trees.
    int posX, posY;

    // Used as storage for computing median, diameter, etc.
    float orScore = 0;   // This event appears in some tree
    float andScore = 0;  // This event appears in both trees
  };

  /// Mapping node class for reconciliation graphs.
  /// Essentially just a pair of one parasite vertex and one host vertex,
  /// along with any children event nodes.
  class MappingNode {
  public:
    Vertex* paraVertex;
    Vertex* hostVertex;
    list<EventNode*> children;
    list<EventNode*> parents;

    MappingNode(Vertex* p, Vertex* h)
    : paraVertex(p), hostVertex(h) { }

    void addChild(EventNode* e) {
        children.push_back(e);
    }

    // true iff this node actually appears in the reconciliaton graph
    bool hit = false; 

    // Utilities
    int dpScore;
    double bottomUpScore; // Used in freqeuncy calulation
    double numberOfTrees; // Number or trees appeared in
    float frequency;

    string str() {
        return "(" + paraVertex->label + "," + hostVertex->label + ")";
    }

    string cstr() { // Condensed string
        return "" + paraVertex->label + "" + hostVertex->label + "";
    }

    // Used as storage for drawing trees.
    int posX, posY;

    // Used to help with loss-reachable calculations
    EventNode* leftLossEvent;
    EventNode* rightLossEvent;
  };
};
typedef Node::EventNode EventNode; // Garbage defs so that it compiles...
typedef Node::MappingNode MappingNode;

/// Important! Here are the representations of reconciliatiton trees/graphs
typedef vector<vector<MappingNode*>> ReconciliationGraph;
typedef list<EventNode*> ReconciliationTree;

//-------|---------|---------|---------|---------|---------|---------|---------|
/// This class stores all of the information you need for solving a DTL
/// problem.
class DTLInstance {
public:
    Tree* hostTree;
    Tree* paraTree;
    unordered_map<Vertex*, Vertex*> tipMapping;
    size_t costD, costT, costL, costS;

    DTLInstance() = default;
    DTLInstance(const DTLInstance&) = default;

    DTLInstance(Tree* ht, Tree* pt, unordered_map<Vertex*, Vertex*> tm) 
        : hostTree(ht), paraTree(pt), tipMapping(tm),
          costD(0), costT(0), costL(0), costS(0) { }

    ~DTLInstance() {
        if (hostTree == nullptr || paraTree == nullptr || mappingTable.size() <= 0) 
            return; // Was probably default constructed
        for (Vertex* vp: paraTree->preOrderVertexList) {
            for (Vertex* vh: hostTree->preOrderVertexList) {
                MappingNode* m = mappingTable[vp->index][vh->index];
                for (EventNode* e: m->children)
                    if (e) delete e;
                if (m) delete m;
            }
        }
        if (hostTree) delete hostTree;
        if (paraTree) delete paraTree;  
    }

    /// Important. Used to set the cost of certain events
    void setCosts(size_t D, size_t T, size_t L, size_t S) {
        costD = D;
        costT = T;
        costL = L;
        costS = S;
    }

    void printInfo() {
        cout << "Here is a DTL problem... Host size: " 
             << SSTR(hostTree->size) + " Para size: " 
             << SSTR(paraTree->size) << endl;
        cout << hostTree->printToString() << endl;
        cout << paraTree->printToString() << endl;
        for (Vertex* v: paraTree->vertexList) {
            if (v->isLeaf()) {
                cout << v->label + ":" + tipMapping[v]->label << ",";
            }
        }
        cout << endl;
    }

    /// Does everything and gives a nice comma delimited string
    void computeAllStatistics(int D, int T, int L, string& csvOut, bool print) {
        setCosts(D, T, L, 0);
        dp();                           if (print) cout << "dp done...";
        createReconciliationGraph();    if (print) cout << "rec graph created...";
        computeNodeFrequency();         if (print) cout << "freq done...";
        computeOtherStatistics();       if (print) cout << "stats done...";

        float asymMed, symMed, symRadius;
        rtAM = computeAsymmetricMedian(asymMed);
        rtSM = computeSymmetricMedian(symMed);
        rtFSM = findFurthestTree(rtSM, symRadius);
        if (print) cout << "medians done...";

        float diameter = computeDiameter();                    if (print) cout << "diameter done...";
        float diameterLossLess = computeDiameterLossLess();    if (print) cout << "lossless diameter done...";
        float minIntersection = computeMinimumIntersection();   if (print) cout << "minintersection done...";
    
        if (print) cout << "DONE!!" << endl;

        csvOut = "";
        csvOut += "," + SSTR(D);
        csvOut += "," + SSTR(T);
        csvOut += "," + SSTR(L);
        csvOut += "," + SSTR(hostTreeSize);
        csvOut += "," + SSTR(hostTreeTipCount);
        csvOut += "," + SSTR(paraTreeSize);
        csvOut += "," + SSTR(paraTreeTipCount);
        csvOut += "," + SSTR(optimalScore);
        csvOut += "," + SSTR(reconciliaitonTreeCount);
        csvOut += "," + SSTR(mappingNodeCount);

        csvOut += "," + SSTR(eventNodeCount);
        csvOut += "," + SSTR(cEventNodeCount);
        csvOut += "," + SSTR(dEventNodeCount);
        csvOut += "," + SSTR(tEventNodeCount);
        csvOut += "," + SSTR(lEventNodeCount);
        csvOut += "," + SSTR(sEventNodeCount);
        csvOut += "," + SSTR(avgEventNodesPerMappingNode);
        csvOut += "," + SSTR(maxEventNodesPerMappingNode);

        csvOut += "," + SSTR(diameter);
        csvOut += "," + SSTR(diameterLossLess);
        csvOut += "," + SSTR(minIntersection);
        csvOut += "," + SSTR(numberOfFreq1Events);
        csvOut += "," + SSTR(numberOfFreq1DTSEvents);

        csvOut += "," + SSTR(asymMed);
        csvOut += "," + SSTR(symMed);
        csvOut += "," + SSTR(symRadius);

        if (print) cout << csvOut << endl;
    }

    // DP stuff
    int INF = 1000000000;
    int optimalScore;

    // This is a vector<vector<MappingNode*>>. It contains all
    // possible mapping nodes in post order for each tree. Some of these
    // mapping nodes may not be used. Run createReconciliationGraph() first,
    // then use m->hit to see if a mapping node is actually hit.
    // The inner vector represents each "group" in the tree.
    ReconciliationGraph mappingTable;

    // Stores the total number of reconciliation trees
    // It's a double since it can exceed 10^18
    double reconciliaitonTreeCount;

    // True if our stuff cannot fit inside of a double.
    // Anything that uses frequency is not accurate anymore.
    bool tooManyReconciliationTrees = false;

    /// This function runs the DP algorithm.
    /// Currently, it uses the semi-optimized O(n^2(logn + 1 + t)) implementation.
    /// That is, for each mapping node, we will look at only the best possible
    /// transfers, then jump out when the score is too high. To do that, we
    /// presort ransfers by score. 
    void dp() {

        optimalScore = INF;

        for (Vertex* vp: paraTree->vertexList) { // post order

            vector<MappingNode*> thisGroup = vector<MappingNode*>(); 

            // Presort the transfers, it's a lot faster.
            // Technically hurts run-time compared to the best-switch
            // stuff. But this was easier to implement.
            // Also, it is still probably faster because you don't
            // need a map.
            vector<MappingNode*> leftTransfers, rightTransfers;
            if (!vp->isLeaf()) {
                leftTransfers = mappingTable[vp->left->index];
                rightTransfers = mappingTable[vp->right->index];
                sort(leftTransfers.begin(), leftTransfers.end(), 
                    [ ]( MappingNode*& lhs, MappingNode*& rhs ) {
                        return lhs->dpScore < rhs->dpScore;
                    });
                sort(rightTransfers.begin(), rightTransfers.end(), 
                    [ ]( MappingNode*& lhs, MappingNode*& rhs ) {
                        return lhs->dpScore < rhs->dpScore;
                    });
            }

            for (Vertex* vh: hostTree->vertexList) { // post order
                
                int bestScore = INF;
                list<EventNode> candidate = list<EventNode>();
                MappingNode* m = new MappingNode(vp, vh);
                thisGroup.push_back(m);
                
                // C event... Just check if the tipMapping is OK
                if (vp->isLeaf() && tipMapping[vp] == vh) {
                    EventNode e(EventType::C, nullptr, nullptr, m);
                    e.dpScore = 0;

                    bestScore = min<int>(bestScore, e.dpScore);
                    candidate.push_back(e);
                }

                // L event... Stay on this group
                if (!vh->isLeaf() && !vp->isRoot()) {
                    // left loss... move down left on the host tree
                    MappingNode* l1 = thisGroup[vh->left->index];
                    EventNode e1(EventType::L, l1, nullptr, m);
                    e1.dpScore = l1->dpScore + costL;

                    // right loss... move down right on the host tree
                    MappingNode* l2 = thisGroup[vh->right->index];
                    EventNode e2(EventType::L, nullptr, l2, m);
                    e2.dpScore = l2->dpScore + costL;

                    bestScore = min<int>(bestScore, e1.dpScore);
                    candidate.push_back(e1);
                    bestScore = min<int>(bestScore, e2.dpScore);
                    candidate.push_back(e2);
                }

                // D event... Don't move on the host tree, but branch on para
                if (!vp->isLeaf()) {
                    MappingNode* d1 = mappingTable[vp->left->index][vh->index];
                    MappingNode* d2 = mappingTable[vp->right->index][vh->index];
                    EventNode e(EventType::D, d1, d2, m);
                    e.dpScore = d1->dpScore + d2->dpScore + costD;

                    bestScore = min<int>(bestScore, e.dpScore);
                    candidate.push_back(e);
                }

                // S event... Branch along both ways
                if (!vp->isLeaf() && !vh->isLeaf()) {
                    // cis speciation
                    MappingNode* s1 = mappingTable[vp->left->index][vh->left->index];
                    MappingNode* t1 = mappingTable[vp->right->index][vh->right->index];
                    EventNode e1(EventType::S, s1, t1, m);
                    e1.dpScore = s1->dpScore + t1->dpScore + costS;

                    // trans speciation
                    MappingNode* s2 = mappingTable[vp->left->index][vh->right->index];
                    MappingNode* t2 = mappingTable[vp->right->index][vh->left->index];
                    EventNode e2(EventType::S, s2, t2, m);
                    e2.dpScore = s2->dpScore + t2->dpScore + costS;

                    bestScore = min<int>(bestScore, e1.dpScore);
                    candidate.push_back(e1);
                    bestScore = min<int>(bestScore, e2.dpScore);
                    candidate.push_back(e2);
                }

                // T event... One stays, the other goes anywhere it can
                if (!vp->isLeaf()) {
                    // left stays
                    MappingNode* d1 = mappingTable[vp->left->index][vh->index];
                    for (MappingNode* t1: rightTransfers) {
                        Vertex* vt = t1->hostVertex;
                        if (hostTree->ancestorTable[vh->index][vt->index] != Rel::INCOMPARABLE)
                            continue; // Illegal transfer...
                        EventNode e1(EventType::T, d1, t1, m);
                        e1.dpScore = d1->dpScore + t1->dpScore + costT;
                        if (e1.dpScore > bestScore) break; // the transfer events
                        // are sorted, so we can just kick out here, we won't
                        // see anything better than this.

                        bestScore = min<int>(bestScore, e1.dpScore);
                        candidate.push_back(e1);
                    }
                    
                    // right stays
                    MappingNode* d2 = mappingTable[vp->right->index][vh->index];
                    for (MappingNode* t1: leftTransfers) {
                        Vertex* vt = t1->hostVertex;
                        if (hostTree->ancestorTable[vh->index][vt->index] != Rel::INCOMPARABLE)
                            continue; // Illegal transfer...
                        EventNode e1(EventType::T, t1, d2, m);
                        e1.dpScore = d2->dpScore + t1->dpScore + costT;
                        if (e1.dpScore > bestScore) break;

                        bestScore = min<int>(bestScore, e1.dpScore);
                        candidate.push_back(e1);
                    }
                }

                // Collect best events from all candidates;
                // they become the children of our mapping node
                for (EventNode e: candidate) {
                    if (e.dpScore == bestScore && e.dpScore < INF) {
                        m->addChild(new EventNode(e.type, e.left, e.right, e.parent));
                    }
                }
                m->dpScore = bestScore;

                if (vp->isRoot()) {
                    optimalScore = min<int>(bestScore, optimalScore);
                }

            }

            mappingTable.push_back(thisGroup);
        }

        // "hit" optimal root mapping nodes so they are marked for later
        for (MappingNode* m: mappingTable[paraTree->root->index]) {
            if (m->dpScore == optimalScore) {
                m->hit = true;
            }
        }
    }

    /// Generate a reconciliation graph from the mapping table
    /// This means throwing away all of the nodes that don't get hit
    /// Use a pre order traversal to see what gets hit.
    /// runtime O(GS(1+t))
    void createReconciliationGraph() {
        for (Vertex* vp: paraTree->preOrderVertexList) {
            for (Vertex* vh: hostTree->preOrderVertexList) {
                MappingNode* m = mappingTable[vp->index][vh->index];
                m->leftLossEvent = nullptr;
                m->rightLossEvent = nullptr;
                if (m->hit) {
                    for (EventNode* e: m->children) {
                        if (e->left != nullptr) {
                            e->left->hit = true;
                            e->left->parents.push_back(e);
                        }
                        if (e->right != nullptr) {
                            e->right->hit = true;
                            e->right->parents.push_back(e);
                        }
                        if (e->type == EventType::L) {
                            if (e->left != nullptr)
                                m->leftLossEvent = e;
                            if (e->right != nullptr)
                                m->rightLossEvent = e;
                        }
                    }
                } else {
                    for (EventNode* e: m->children) {
                        if (e) delete e; // get rid of non-hit event nodes
                    }
                    m->children.clear();
                }
            }
        }
    }

    /// Computes the frequency of each mapping and event node. 
    /// I.e. the number of reconciliation trees that a node appears in.
    /// This function also counts the number of total reconciliation trees.
    /// runtime: O(GS(1+t))
    void computeNodeFrequency() {
        reconciliaitonTreeCount = 0;
        // This is the bottom up step, assigning a bottom up score to each
        // node. That represents the total number of reconciliation subtrees
        // rooted at that node.
        for (Vertex* vp: paraTree->vertexList) { // post order
            for (Vertex* vh: hostTree->vertexList) { // post order
                MappingNode* m = mappingTable[vp->index][vh->index];
                if (!m->hit) continue;

                double sum = 0;
                for (EventNode* e: m->children) {
                    double prod = 1;
                    if (e->left != nullptr)
                        prod *= e->left->bottomUpScore;
                    if (e->right != nullptr)
                        prod *= e->right->bottomUpScore;
                    e->bottomUpScore = prod;
                    e->numberOfTrees = 0;
                    e->frequency = 0;
                    sum += e->bottomUpScore;
                }
                m->bottomUpScore = sum;
                m->numberOfTrees = 0;
                m->frequency = 0;

                if (m->bottomUpScore > 1e150 ) {
                    cout << "Warning: too many reconciliaiton trees." << endl;
                    cout << "There are more than 10^150. " << endl;
                    tooManyReconciliationTrees = true;
                    return;
                }

                // Sum the scores of the mapping nodes that are from the
                // root of the parasite tree
                if (vp->isRoot()) {
                    reconciliaitonTreeCount += m->bottomUpScore;
                }   
            }
        }
        // This is the top down step. Here, the nodes pass their scores down
        // to their kids.
        for (Vertex* vp: paraTree->preOrderVertexList) {
            for (Vertex* vh: hostTree->preOrderVertexList) {
                MappingNode* m = mappingTable[vp->index][vh->index];
                if (!m->hit) continue;
                
                // Root mapping nodes already know their info
                if (vp->isRoot()) {
                    m->numberOfTrees = m->bottomUpScore;
                }

                m->frequency = (float)(m->numberOfTrees / reconciliaitonTreeCount);
                
                for (EventNode* e: m->children) {
                    e->numberOfTrees = (m->numberOfTrees / m->bottomUpScore) * e->bottomUpScore;
                    e->frequency = (float)(e->numberOfTrees / reconciliaitonTreeCount);
                    if (e->left != nullptr)
                        e->left->numberOfTrees += e->numberOfTrees;
                    if (e->right != nullptr)
                        e->right->numberOfTrees += e->numberOfTrees;
                }
            }
        }
    }

    // Other statistics
    int hostTreeSize, paraTreeSize;
    int mappingNodeCount, internalMappingNodeCount, leafMappingNodeCount;
    int eventNodeCount;
    int cEventNodeCount, dEventNodeCount, tEventNodeCount;
    int lEventNodeCount, sEventNodeCount;
    int hostTreeTipCount, paraTreeTipCount;
    float avgEventNodesPerMappingNode;
    int maxEventNodesPerMappingNode;
    int numberOfFreq1Events, numberOfFreq1DTSEvents;

    ReconciliationTree rtAM, rtSM, rtFSM;

    /// computes a number of statistics on the reconciliation graph
    void computeOtherStatistics() {
    
        hostTreeSize = hostTree->size; 
        paraTreeSize = paraTree->size;

        mappingNodeCount = 0;
        internalMappingNodeCount = 0;
        leafMappingNodeCount = 0;
        eventNodeCount = 0;
        cEventNodeCount = 0;
        dEventNodeCount = 0;
        tEventNodeCount = 0;
        lEventNodeCount = 0;
        sEventNodeCount = 0;
        maxEventNodesPerMappingNode = 0;
        numberOfFreq1DTSEvents = 0;
        numberOfFreq1Events = 0;

        for (Vertex* vp: paraTree->vertexList) { // post order
            for (Vertex* vh: hostTree->vertexList) { // post order
                MappingNode* m = mappingTable[vp->index][vh->index];
                if (!m->hit) continue;
                mappingNodeCount++;
                if (vp->isLeaf()) 
                    leafMappingNodeCount++;
                else internalMappingNodeCount++;

                for (EventNode* e: m->children) {
                    if (e->frequency >= 1) {
                        numberOfFreq1Events++;
                        if (e->type != EventType::L && e->type != EventType::C)
                            numberOfFreq1DTSEvents++;
                    }
                    eventNodeCount++;
                    switch(e->type) {
                    case EventType::C: cEventNodeCount++; break; 
                    case EventType::D: dEventNodeCount++; break; 
                    case EventType::T: tEventNodeCount++; break; 
                    case EventType::L: lEventNodeCount++; break; 
                    case EventType::S: sEventNodeCount++; break;
                    }
                }
                maxEventNodesPerMappingNode = max<int>(m->children.size(),
                    maxEventNodesPerMappingNode);
            }
        }

        hostTreeTipCount = (hostTreeSize+1)/2;
        paraTreeTipCount = (paraTreeSize+1)/2;
        avgEventNodesPerMappingNode = (eventNodeCount * 1.0f) / mappingNodeCount;
    }

    /// For sorting
    bool sortByDpScore(MappingNode*& m1, MappingNode*& m2) {
        return m1->dpScore < m2->dpScore;
    }

    /// Print out the nodes in our reconciliation graph
    void printReconciliationGraph() {
        cout << "Here comes a reconciliation graph... " << endl;
        for (Vertex* vp: paraTree->vertexList) { // post order
            for (Vertex* vh: hostTree->vertexList) { // post order
                MappingNode* m = mappingTable[vp->index][vh->index];
                if (!m->hit) continue; // Ignore not hit nodes
                cout << "Mapping Node: " << m->str() 
                    << " Score: " << SSTR(m->dpScore) 
                    << " Freq: " << SSTR(m->frequency) << endl;
                for (EventNode* e: m->children)       
                    cout << "---Event: " << e->str()
                         << " Freq: " << SSTR(e->frequency) << endl;
            }
        }
        cout << "This reconciliation graph has #trees = " 
             << reconciliaitonTreeCount << endl;
    }

    /// Print a reconciliaiton tree
    void printReconciliationTree(ReconciliationTree rt) {
        cout << "Here is a reconciliation tree... " << endl;
        for (EventNode* e: rt) {
            cout << e->str() << endl;
        }
    }

    //----------------------ONE TREE--------------------------------------

    /// This function will find the reconciliation tree R that maximizes the sum
    /// over all event nodes E in R of E.orScore. Thus, before you call this
    /// function, you should set orScore to the way you want each event to be
    /// scored. For example, use frequency to get the asymmetric median.
    /// Use negative scores to minimize something.
    /// It returns both the max score, and one arbitrary R-tree that attains it.
    ReconciliationTree maximizeScoreOneTree(/*Out*/ float& maximumScore) {
        // Note! orScore should hold the desired values by this point!

        unordered_map<MappingNode*, float> table;
        unordered_map<MappingNode*, EventNode*> rTable; // for reconstruction

        for (Vertex* vp: paraTree->vertexList) { // post order
            for (Vertex* vh: hostTree->vertexList) { // post order
                MappingNode* m = mappingTable[vp->index][vh->index];
                if (!m->hit) continue; // Ignore not hit nodes
                float bestScore = -INF;
                EventNode* bestEvent;

                for (EventNode* e: m->children) {
                    float thisScore = e->orScore;
                    if (e->left != nullptr) 
                        thisScore += table[e->left];
                    if (e->right != nullptr) 
                        thisScore += table[e->right];
                    if (thisScore > bestScore) {
                        bestScore = thisScore;
                        bestEvent = e;
                    }
                }

                table[m] = bestScore;
                rTable[m] = bestEvent;
            }
        }

        // loop over all root mapping nodes to find the best
        float bestScore = -INF;
        MappingNode* bestMappingNode;
        for (Vertex* vh: hostTree->preOrderVertexList) {
            MappingNode* m = mappingTable[paraTree->root->index][vh->index];
            if (!m->hit) continue; // Ignore not hit nodes
            
            if (table[m] > bestScore) {
                bestScore = table[m];
                bestMappingNode = m;
            }
        }
        maximumScore = bestScore;

        // Reconstruct reconciliation tree
        list<EventNode*> rTree;
        queue<MappingNode*> mQueue;
        mQueue.push(bestMappingNode);
        while (!mQueue.empty()) {
            MappingNode* m = mQueue.front();
            mQueue.pop();
            EventNode* e = rTable[m];
            rTree.push_back(e);
            if (e->left != nullptr) 
                mQueue.push(e->left);
            if (e->right != nullptr) 
                mQueue.push(e->right);
        }

        return rTree;
    }

    /// Returns the asymmetric median
    ReconciliationTree computeAsymmetricMedian(/*Out*/ float& maximumFreq) {
        for (Vertex* vp: paraTree->vertexList) {
            for (Vertex* vh: hostTree->vertexList) {
                MappingNode* m = mappingTable[vp->index][vh->index];
                for (EventNode* e: m->children) {
                    e->orScore = e->frequency;
                }
            }
        }

        return maximizeScoreOneTree(maximumFreq);
    }

    /// Returns the symmetric median
    ReconciliationTree computeSymmetricMedian(/*Out*/ float& maximumFreq) {
        for (Vertex* vp: paraTree->vertexList) {
            for (Vertex* vh: hostTree->vertexList) {
                MappingNode* m = mappingTable[vp->index][vh->index];
                for (EventNode* e: m->children) {
                    e->orScore = e->frequency - 0.5f;
                }
            }
        }

        return maximizeScoreOneTree(maximumFreq);
    }

    /// Finds the reconciliation tree that is furthest away by symmetric
    /// distance from the given reconciliation tree.
    ReconciliationTree findFurthestTree(ReconciliationTree rt, /*Out*/ float& maximumDist) {
        for (Vertex* vp: paraTree->vertexList) {
            for (Vertex* vh: hostTree->vertexList) {
                MappingNode* m = mappingTable[vp->index][vh->index];
                for (EventNode* e: m->children) {
                    e->orScore = 1;
                }
            }
        }
        int nodeCount = 0;
        for (EventNode* e: rt) {
            e->orScore = -1;
            nodeCount++;
        }

        ReconciliationTree ret = maximizeScoreOneTree(maximumDist);
        maximumDist += nodeCount;
        return ret;
    }

    //----------------------TWO TREES-------------------------------------

    /// This function finds the pair of reconcilation trees that maximizes
    /// the orScore over all event nodes contained in either tree plus
    /// the andScore over all event nodes contained in both trees. These values
    /// need to be set before this function is called. 
    float maximizeScoreTwoTrees() {

        // This is the enter mapping table. It is a 3 dimensional table.
        // The first entry is the group. The second and third entry are the 
        // pair of mapping nodes considered.
        // It stores the best score for any two reconciliation sub-trees
        // rooted at those two mapping nodes.
        vector<vector<vector<float>>> enter (paraTree->size, 
                vector<vector<float>>(hostTree->size, 
                vector<float>(hostTree->size, -INF)));
        // Note: the hacky vector and index stuff is faster than a
        // unordered_map by quite a bit
        
        // This is the playnice mapping table. It is a 3 dimensional table.
        // The first entry is the group. The second and third entry are the 
        // pair of mapping nodes considered.
        // It stores the best score for any two reconciliation sub-trees
        // rooted at those two mapping nodes, subject to the additional
        // constraint that one of the mapping nodes is a direct ancestor of
        // the other and the ancestor "plays nice" by immediately exiting
        // the group so that the other node is free to take loss events. 
        vector<vector<vector<float>>> playnice (paraTree->size, 
                vector<vector<float>>(hostTree->size, 
                vector<float>(hostTree->size, -INF)));

        for (Vertex* vp: paraTree->vertexList) { // post order
            //if (vp->index%100==0) cout << SSTR(vp->index) << endl;
            for (Vertex* vh1: hostTree->vertexList) { // post order
                for (Vertex* vh2: hostTree->vertexList) { // post order

                    MappingNode* m1 = mappingTable[vp->index][vh1->index];
                    MappingNode* m2 = mappingTable[vp->index][vh2->index];
                    if (!m1->hit || !m2->hit) continue;

                    // What we will be computing here. Initially -INF
                    float& thisEntry = enter[vp->index][vh1->index][vh2->index];

                    // In any situation, we can take exit events from both
                    // m1 and m2. So, let's do that first.
                    for (EventNode* e1: m1->children) {
                        for (EventNode* e2: m2->children) {
                            if (e1->type == EventType::L || e2->type == EventType::L)
                                continue; // Not an exit event

                            // Trivial C event, occurs only if vp is a leaf
                            if (e1->type == EventType::C && e1 == e2)
                                thisEntry = e1->andScore + e1->orScore;

                            // Otherwise, we have two exit-events (D, T, or S)
                            else {
                                float eq = e1 == e2 ? e1->orScore + e1->andScore : e1->orScore + e2->orScore;
                                float cL = enter[vp->left->index][e1->left->hostVertex->index][e2->left->hostVertex->index];
                                float cR = enter[vp->right->index][e1->right->hostVertex->index][e2->right->hostVertex->index];
                                thisEntry = max<float>(eq+cL+cR, thisEntry);
                            }
                        }
                    }
                    
                    // m1 and m2 are incomparable
                    if (hostTree->ancestorTable[vh1->index][vh2->index] == Rel::INCOMPARABLE ) {
                        
                        if (m1->leftLossEvent != nullptr) {
                            EventNode* e = m1->leftLossEvent;
                            float ls = e->orScore + enter[vp->index][e->left->hostVertex->index][vh2->index];
                            thisEntry = max<float>(thisEntry, ls);
                        }
                        if (m1->rightLossEvent != nullptr) {
                            EventNode* e = m1->rightLossEvent;
                            float ls = e->orScore + enter[vp->index][e->right->hostVertex->index][vh2->index];
                            thisEntry = max<float>(thisEntry, ls);
                        }
                        if (m2->leftLossEvent != nullptr) {
                            EventNode* e = m2->leftLossEvent;
                            float ls = e->orScore + enter[vp->index][vh1->index][e->left->hostVertex->index];
                            thisEntry = max<float>(thisEntry, ls);
                        }
                        if (m2->rightLossEvent != nullptr) {
                            EventNode* e = m2->rightLossEvent;
                            float ls = e->orScore + enter[vp->index][vh1->index][e->right->hostVertex->index];
                            thisEntry = max<float>(thisEntry, ls);
                        }
                    }

                    // m1 and m2 are equal
                    else if (m1 == m2) {

                        if (m1->leftLossEvent != nullptr) {
                            EventNode* e = m1->leftLossEvent;
                            float ls = e->andScore + e->orScore
                                    + enter[vp->index][e->left->hostVertex->index][e->left->hostVertex->index];
                            thisEntry = max<float>(thisEntry, ls);

                            float lp = e->orScore
                                    + playnice[vp->index][vh1->index][e->left->hostVertex->index];
                            thisEntry = max<float>(thisEntry, lp);
                        }
                        if (m1->rightLossEvent != nullptr) {
                            EventNode* e = m1->rightLossEvent;
                            float ls = e->andScore + e->orScore
                                    + enter[vp->index][e->right->hostVertex->index][e->right->hostVertex->index];
                            thisEntry = max<float>(thisEntry, ls);

                            float lp = e->orScore
                                    + playnice[vp->index][vh1->index][e->right->hostVertex->index];
                            thisEntry = max<float>(thisEntry, lp);
                        }
                        if (m1->leftLossEvent != nullptr && m1->rightLossEvent != nullptr) {
                            EventNode* el = m1->leftLossEvent;
                            EventNode* er = m1->rightLossEvent;
                            float ls = el->orScore + er->orScore
                                    + enter[vp->index][el->left->hostVertex->index][er->right->hostVertex->index];
                            thisEntry = max<float>(thisEntry, ls);
                        }
                    }

                    // m1 and m2 are comparable
                    else {
                        float& thisPlaynice = playnice[vp->index][vh1->index][vh2->index];
                        thisPlaynice = max<float>(thisPlaynice, thisEntry);

                        if (hostTree->ancestorTable[vh1->index][vh2->index] == Rel::ANCESTOR ) {
                            // cheekily steal the symmetric value
                            thisPlaynice = playnice[vp->index][vh2->index][vh1->index];
                            thisEntry = enter[vp->index][vh2->index][vh1->index];
                            continue;
                        } 
                        
                        // At this point, vh2 is the ancestor of vh1
                        if (m1->leftLossEvent != nullptr) {
                            EventNode* e = m1->leftLossEvent;
                            float ls = e->orScore
                                    + playnice[vp->index][e->left->hostVertex->index][vh2->index];
                            thisPlaynice = max<float>(thisPlaynice, ls);
                        }
                        if (m1->rightLossEvent != nullptr) {
                            EventNode* e = m1->rightLossEvent;
                            float ls = e->orScore
                                    + playnice[vp->index][e->right->hostVertex->index][vh2->index];
                            thisPlaynice = max<float>(thisPlaynice, ls);
                        }

                        thisEntry = max<float>(thisEntry, thisPlaynice);
                        if (m2->leftLossEvent != nullptr) {
                            EventNode* e = m2->leftLossEvent;
                            float ls = e->orScore + enter[vp->index][vh1->index][e->left->hostVertex->index];
                            thisEntry = max<float>(thisEntry, ls);
                        }
                        if (m2->rightLossEvent != nullptr) {
                            EventNode* e = m2->rightLossEvent;
                            float ls = e->orScore + enter[vp->index][vh1->index][e->right->hostVertex->index];
                            thisEntry = max<float>(thisEntry, ls);
                        }
                    }
                }
            }
        }

        // Find the best ones among all root pairs
        float best = -INF;
        for (Vertex* vh1: hostTree->vertexList) {
            for (Vertex* vh2: hostTree->vertexList) {
                best = max<float>(best, 
                    enter[paraTree->root->index][vh1->index][vh2->index]);
            }
        }
        return best;
    }


    /// Returns the diameter, the distance between the two furthest apart trees
    /// using symmetric set scoring
    float computeDiameter() {
        for (Vertex* vp: paraTree->vertexList) {
            for (Vertex* vh: hostTree->vertexList) {
                MappingNode* m = mappingTable[vp->index][vh->index];
                for (EventNode* e: m->children) {
                    e->orScore = 1;
                    e->andScore = -1;
                }
            }
        }

        return maximizeScoreTwoTrees();
    }

    /// Returns the diameter, the distance between the two furthest apart trees
    /// using symmetric set scoring, not counting loss events
    float computeDiameterLossLess() {
        for (Vertex* vp: paraTree->vertexList) {
            for (Vertex* vh: hostTree->vertexList) {
                MappingNode* m = mappingTable[vp->index][vh->index];
                for (EventNode* e: m->children) {
                    e->orScore = 1;
                    e->andScore = -1;
                    if (e->type == EventType::L) {
                        e->orScore = 0;
                        e->andScore = 0;
                    }
                }
            }
        }

        return maximizeScoreTwoTrees();
    }

    /// Returns the diameter, the distance between the two furthest apart trees
    /// using symmetric set scoring
    float computeMinimumIntersection() {
        for (Vertex* vp: paraTree->vertexList) {
            for (Vertex* vh: hostTree->vertexList) {
                MappingNode* m = mappingTable[vp->index][vh->index];
                for (EventNode* e: m->children) {
                    e->orScore = 0;
                    if (e->type != EventType::C)
                        e->andScore = -1;
                    else e->andScore = 0;
                }
            }
        }

        return -1.0f * maximizeScoreTwoTrees();
    }


    //----------------------K COVER-------------------------------------

    // The approximation algorithm described in the 2015 DTLRnB paper
    // Take the best remaining tree, then set the score of all the events
    // found to zero to not see them again. Approximates the k best trees
    // that maximize collective frequency.
    // Runtime: O(kGSt)
    float kCoverWMCP(size_t k, bool includeLoss, /*Out*/ list<ReconciliationTree>& rts) {
        for (Vertex* vp: paraTree->vertexList) {
            for (Vertex* vh: hostTree->vertexList) {
                MappingNode* m = mappingTable[vp->index][vh->index];
                for (EventNode* e: m->children) {
                    if (e->type == EventType::L && !includeLoss)
                        e->orScore = 0;
                    else e->orScore = e->frequency;
                }
            }
        }
        float score = 0;

        for (size_t i = 0; i < k; i++) {
            float thisScore;
            ReconciliationTree rt = maximizeScoreOneTree(thisScore);
            score += thisScore;
            cout << SSTR(score) << endl;
            rts.push_back(rt);
            for (EventNode* e: rt) {
                e->orScore = 0;
            }
        }

        return score;
    }


};

//-------|---------|---------|---------|---------|---------|---------|---------|
/// A class for parsing stuff
class Parser {
public:

    /// This function parses a newick tree
    Tree* parseTree(string s) {
        pair<Vertex*, string> res = parseTreeHelper(s);
        Tree* t = new Tree(res.first);
        return t;
    }

    /// This function takes a file name and parses
    /// two trees and a tip mapping into a DTL problem
    DTLInstance parseFile(string fileName) {
        ifstream file;
        file.open(fileName);
        if(!file.is_open()) {
            throw invalid_argument("Bad File");
        }
    
        // parse host and parasite trees
        string fileLine;
        getline(file, fileLine);
        Tree* hostTree = parseTree(fileLine);
        getline(file, fileLine);
        Tree* paraTree = parseTree(fileLine);

        // parse the tip mapping
        unordered_map<Vertex*, Vertex*> tipMapping;
        while(getline(file, fileLine)) {
            parseTipMapping(tipMapping, fileLine, hostTree, paraTree);
        }

        return DTLInstance(hostTree, paraTree, tipMapping);
    }

private:
    /// Here, we reculsively parse a single vertex at a time
    pair<Vertex*, string> parseTreeHelper(string s) {
        if (s.length() == 0) {
            return pair<Vertex*, string>(nullptr, "");
        } else if (s[0] == '(') {
            // If we see a '(', we expect two vertices to appear.
            // Start the substrings at 1 to avoid ( and , chars
            pair<Vertex*, string> lf = parseTreeHelper(s.substr(1));
            pair<Vertex*, string> rt = parseTreeHelper(lf.second.substr(1));
            // Look for the end of the label, which will be ) or ,
            size_t q = 1;
            for (; q < rt.second.length(); q++) {
                if (rt.second[q] == ')' || rt.second[q] == ',') break;
            }
            // Create our new vertex, and we know its children already
            Vertex* v = new Vertex(lf.first, rt.first, rt.second.substr(1,q-1));
            return pair<Vertex*, string>(v, rt.second.substr(q));
        } else {
            size_t q = 0;
            for (; q < s.length(); q++) {
                if (s[q] == ')' || s[q] == ',') break;
            }
            // This is a leaf vertex
            Vertex* v = new Vertex(s.substr(0,q));
            return pair<Vertex*, string>(v, s.substr(q));
        }
    }

    /// Parse a singular tip mapping of the form "a:B" and add it to our map
    void parseTipMapping(unordered_map<Vertex*, Vertex*>& tipMapping, string s,
                            Tree* ht, Tree* pt) {
        size_t colon = s.find(":");
        if (colon == string::npos) return;
        string sp = s.substr(0,colon);
        string sh = s.substr(colon+1);
        Vertex* vp = pt->labelMap[sp];
        Vertex* vh = ht->labelMap[sh];
        tipMapping[vp] = vh;
    }
};

