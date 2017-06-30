### **How to use this code from the command line:**

After typing `python DTLReconGraph.py`, one should enter the values for:

  * The name of the file in which the data are stored and
  * The costs of duplication, transfer, and loss, 
  
in that order. Optionally, one could include (anywhere in the arguments) `-s` for the code to use a slower computation method but which would then result in the code including frequency scoring of events in the output (see explanation of the output for more on this). If the user doesn't give input in the expected form, a usage statement will be displayed showing the user what inputs to enter and how to enter them. Note that the data file should include: 
  * Two newick trees (formatted, for example, as ((C,D)B,E)A if A has children B and E and B has children C and D) on separate lines, with the species/host tree given first and the gene/parasite tree second, and
  * Mappings between all of the tips of the two trees (formatted, for example, c:C; d:D; e:E if c maps to C, d to D and e to E), with a different one on each line.

This code then prints: 
  * Two Python dictionaries in the format {(parent, child):(parent, child, child 1 of child, child 2 of child):, ...}, where parent is a node in the tree, child is one of parent's children, and children 1 and 2 of child represent children of the child of parent (in other words, grandchildren of the parent) - one for each of the two trees used as input. This is essentially a dictionary representation of the two trees given as input;
  * A Python dictionary with mapping nodes as keys and lists of event nodes that apply to those mapping nodes in a Maximum Parsimony Reconciliation as values, which has the format {(gene node, species node): [(event string, resulting mapping node 1, resulting mapping node 2, optionally frequency scoring), minimum cost], ...}. In this representation, the gene and species node values represent the nodes in the given newick trees that are being mapped to each other in the reconciliation; the event string is a character ('D', 'T', 'L', 'C', or 'S', representing duplication, transfer, loss, contemporary event, and speciation, respectively) representing the event that those mapped nodes undergo; the resulting mapping nodes are the mappings that are induced in the next 'leve' of the reconciliation graph as a result of the given event; lastly, the frequency scoring is the number of times the given event occurs in MPRs for the given data input - this is the extra information that is shown if the user chooses the `-s` option in the command line arguments. Lastly, the minimum cost is the minimum cost associated with this mapping in an MPR;
  * The number of Maximum Parsimony Reconciliations for the given data set and costs, as an integer

Note each of these printed values is separated by a semicolon and two newlines, so the output is easily readable by a user (the newlines) but also contains characters that make the output easily separable, say, by a computer (the semicolons).

### **How to use this code in interactive mode:**

One should simply type `python -i DTLReconGraph.py` to access the contents of the code. From there, all major and helper functions will be available, however the two most important are DP and reconcile. DP is the workhorse of the code - it utilizes several helper functions to actually perform computations and implement the algorithm given in the technical report titled "HMC CS Technical Report CS-2011-1: Faster Dynamic Programming Algorithms for the Cophylogeny Reconstruction Problem". For more details on the algorithm itself, see that report and/or the comments/docstrings included in the file. Given host and parsite trees, tip mappings, costs for duplication, transfer, and loss, and a boolean value indicating whether the user wishes to include frequency scores in the output, DP implements this algorithm and returns the host and parasite trees in dictionary form, the DTL maximum parsimony reconciliation graph, and the number of reconciliations for the given trees and tip mappings (see the outputs discussed in the section on using this code from the command line for more detail on the output format). reconcile, however, is the more practically useful function. Since the data are implemented as newick trees and mappings, reconcile utilizes a separate module that both handles getting the data from a separate file and reformats the inputs to work nicely with DP. Reconcile reformats the species and gene trees to match the output format given in the section on running from the command line, and prints out these trees along with the reconciliation graph (again, in the format discussed above) and the number of Maximum Parsimony Reconciliations as an integer. Although one may play with any/all functions included in this file and those on which it depends, DP and reconcile are the most important functions.

## __**How to use Diameter.py**__

### From the Command Line:

#### Single Files

To calculate a single diameter on the command line, use the following command:
> python Diameter.py calc file d t l [logfile]

In this command, `file` is the relative path to the newick file you wish to compute the diameter of, `d`, `t`, and `l` are the event costs, and `[logfile]` is the path to the csv file you wish to log to (do not include the `.csv` extension; it will be created for you). If logfile is not specified, no csv will be created.

This argument will calculate the diameter twice, once counting losses as part of the diameter, and once not counting losses. Both results will be printed to the screen, and each result will be recorded to a seperate csv file (the zero loss diameter will have `_zl` at the end of its filename).

#### Many Files

For repeated diameter computations, use the following command:
> python Diameter.py rep files start end d t l [logfile]

This command works the same as the `calc` command, except that `files` is a file path *pattern* where any `#`'s will be replaced with a number, and all files in the range between the number specifed in the `start` and `end` command will have the diameter calculated.

For example, to calculate the diameters for files `tree1.newick`, `tree2.newick`, and `tree4.newick` (all located in a folder called `data`) with DTL costs of 2, 3, and 4, use this command:
>python Diameter.py rep "data/tree#.newick" 1 5 2 3 4

The program will skip over the non-existant `tree3.newick`.

#### Misc

If you like looking at pretty tables, use the ```python Diameter.py test``` command. It will print out some pretty tables.

### Via Interactive Mode:

#### Single Files

To calculate the diameter for a single file in interactive mode, call the following function:
>calculate_diameter_from_file(filename, D, T, L, csv_file=None, debug=False)

The arguments in this function are similar to the command line option `calc`. `filename` is the path to the file that will be reconciled and the DTL Diameter calculated for. `D`, `T`, and `L` are the (integer) event costs, and `csv_file` is the a string containing path to the csv file used for logging (use `None` for no logging).

`debug`, when set to true, will make Diameter print out the dynamic programming tables it uses to find the solution (but only the ones that are less than 30x30 size).

#### Many Files

To calculate a set of numbered files, use this function:
> repeatedly_calculate_diameter(file_pattern, start, end, d, t, l, log=None, debug=False)

Where `file_pattern` is the pattern used to find the right files (as described in the command line section), `start` is the starting file number, `end` is the exclusive ending file number, `d`, `t`, and `l` are the event costs, and `log` is the csv file to log to (or `None` if no logfile is needed, although that is not recommended).

`debug` works the same as in `calculate_diameter_from_file`, but that is *very* not recommended.

