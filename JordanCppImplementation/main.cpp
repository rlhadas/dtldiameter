#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include "tree.cpp"

using namespace std;

int main() {
	cout << "Hello" << endl;

	Parser parser;

	/*DTLInstance dtl = parser.parseFile("TreeLifeData\\atest1.newick");
	dtl.setCosts(0, 0, 0, 0);
    dtl.dp();                          
    dtl.createReconciliationGraph();   
    dtl.computeNodeFrequency();       
    dtl.computeOtherStatistics();
	dtl.printInfo();
	dtl.printReconciliationGraph();

	list<ReconciliationTree> rts;
	cout << "3 Cover LL: " << SSTR(dtl.kCoverWMCP(3, false, rts)) << endl;
	for (ReconciliationTree rt: rts) {
		dtl.printReconciliationTree(rt);
	}

	list<ReconciliationTree> rts2;
	cout << "my 3 Cover LL: " << SSTR(dtl.kCoverMY(3, false, rts)) << endl;
	//for (ReconciliationTree rt: rts2) {
	//	dtl.printReconciliationTree(rt);
	//}
	if (1) return 0;*/

	// .csv line format
	string csv = "File,D,T,L,HostVertexCount,HostTipCount,ParaVertexCount,ParaTipCount,optimalScore,numReconciliationTrees,numMappingNodes,numEventNodes,numCeventNodes,numDeventNodes,numTeventnodes,numLeventnodes,numSeventnodes,avgEventPerMapping,maxEventPerMapping,Diameter,LosslessDiameter,MinIntersection,Freq1Events,Freq1DTLEvent,asymMed,symMed,MaxDistFromMed";
	ofstream myfile;
  	myfile.open("myOutputLog.csv");
  	myfile << csv << endl;
  	
	for (int k = 1; k < 5666; ++k) {
		string intSt = SSTR(k); 
		for (int st = intSt.size(); st<4; ++st) intSt = "0" + intSt;
		cout << intSt << endl;

		try {
			DTLInstance d = parser.parseFile("TreeLifeData\\COG" + intSt + ".newick");
			string s;
			d.computeAllStatistics(1,4,1,s,true);
			myfile << intSt << s << endl;
		} catch (exception& e) {
			// Probably bad file...
			cout << e.what() << endl;
		}
	}

	myfile.close();
	cout << "ALL DONE" << endl;
}
