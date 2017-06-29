# dtldiameter

* __**How to use DTLReconGraph.py**__: Of course, one could run this file in interactive mode and test any functions individually. However, the most important function in this file is DP, especially since its results are used in Diameter.py. The most helpful function to run, however, is reconcile. This is what Diameter.py actually uses to compute results, and it extracts all of the important info from DP and allows the user to pass data files for the program to run on. Consequently, users should focus on understanding the reconcile function and all functions on which it depends in order to get the most out of this program. Of course, all of the other functions that don't directly have anything to do with reconcile would be useful to understand also. All inputs, return values, and their formats (along with the rest of the code) are thoroughly explained in the file.

* __**How to use Diameter.py**__:
