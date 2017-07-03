# dtldiameter

## How to use DTLReconGraph.py

### How to use this code from the command line:

After typing `python DTLReconGraph.py`, one should enter the values for:

  * The name of the file in which the data are stored and
  * The costs of duplication, transfer, and loss, 
  
in that order. Optionally, one could include (anywhere in the arguments) `-s` for the code to use a slower computation method but which would then result in the code including frequency scoring of events in the output (see explanation of the output for more on this). If the user doesn't give input in the expected form, a usage statement will be displayed showing the user what inputs to enter and how to enter them. Note that the data file should include: 
  * Two newick trees (formatted, for example, as ((C,D)B,E)A if A has children B and E and B has children C and D) on separate lines, with the species/host tree given first and the gene/parasite tree second, and
  * Mappings between all of the tips of the two trees (formatted, for example, c:C; d:D; e:E if c maps to C, d to D and e to E), with a different one on each line.

This code then prints: 
  * Two Python dictionaries in the format {(parent, child):(parent, child, child 1 of child, child 2 of child):, ...}, where parent is a node in the tree, child is one of parent's children, and children 1 and 2 of child represent children of the child of parent (in other words, grandchildren of the parent) - one for each of the two trees used as input. This is essentially a dictionary representation of the two trees given as input;
  * A Python dictionary with mapping nodes as keys and lists of event nodes that apply to those mapping nodes in a Maximum Parsimony Reconciliation as values, which has the format {(gene node, species node): [(event string, resulting mapping node 1, resulting mapping node 2, optionally frequency scoring), minimum cost], ...}. In this representation, the gene and species node values represent the nodes in the given newick trees that are being mapped to each other in the reconciliation; the event string is a character ('D', 'T', 'L', 'C', or 'S', representing duplication, transfer, loss, contemporary event, and speciation, respectively) representing the event that those mapped nodes undergo; the resulting mapping nodes are the mappings that are induced in the next 'level' of the reconciliation graph as a result of the given event; lastly, the frequency scoring is the number of times the given event occurs in MPRs for the given data input - this is the extra information that is shown if the user chooses the `-s` option in the command line arguments. Lastly, the minimum cost is the minimum cost associated with this mapping in an MPR;
  * The number of Maximum Parsimony Reconciliations for the given data set and costs, as an integer

Note each of these printed values is separated by a blank line, so the output is easily readable by a user.

### How to use this code in interactive mode:

One should simply type `python -i DTLReconGraph.py` to access the contents of the code. From there, all major and helper functions will be available, however the two most important are DP and reconcile. 

#### DP

DP is the workhorse of the code - it utilizes several helper functions to actually perform computations and implement the algorithm given in the technical report titled "HMC CS Technical Report CS-2011-1: Faster Dynamic Programming Algorithms for the Cophylogeny Reconstruction Problem". For more details on the algorithm itself, see that report and/or the comments/docstrings included in the file. Given host and parasite trees, tip mappings, costs for duplication, transfer, and loss, and a boolean value indicating whether the user wishes to include frequency scores in the output, DP implements this algorithm and returns the host and parasite trees in dictionary form, the DTL maximum parsimony reconciliation graph, and the number of reconciliations for the given trees and tip mappings (see the outputs discussed in the section on using this code from the command line for more detail on the output format).

#### Reconcile

Reconcile is the more practically useful function. Since the data are implemented as newick trees and mappings, reconcile utilizes a separate module that both handles getting the data from a separate file and reformats the inputs to work nicely with DP. Reconcile reformat the species and gene trees to match the output format given in the section on running from the command line, and prints out these trees along with the reconciliation graph (again, in the format discussed above) and the number of Maximum Parsimony Reconciliations as an integer. Although one may play with any/all functions included in this file and those on which it depends, DP and reconcile are the most important functions.

## How to use Diameter.py

### From the Command Line:

#### Required Arguments

Running Diameter from the command line has the following usage pattern:

> Diameter.py [options] file d t l

In this command, `file` is the relative path to the newick file you wish to compute the diameter of, `d`, `t`, and `l` are the event costs. Specifically, `d` is the event cost of duplications, `t` is the event cost of transfers, and `l` is the event cost of losses. This command will calculate the diameter of the file twice, once counting losses as part of the diameter, and once not counting losses. Both results will be printed to the screen.

#### Optional Flags

There are several flags availible:

* `-l` or `--log` takes one argument, and will log the results of the calculation to a row in the provided csv file. Be sure to specify a path without a file extension, as this program will add one for you. Both the regular and zero-loss result will be recorded to seperate files (the zero loss diameter will have `_zl` at the end of its filename).

* `-i` or `--iterate` takes two arguments, and tells the program to iterate over several numbered files and perform the calculation. This first argument is the starting number, and the second argument is the ending number. Note that when this flag is active, the `file` required argument will replace any `#` characters in the file name with whatever number it is currently on, and use that file in its current calculation.

* `-q` or `--quiet` supresses all text output (except that generated by the `-i` flag).

* `-d` or `--debug` prints out every dynamic programming table generated by the Diameter algorithm.

#### Examples

To calculate a single diameter from file `tree1.newick` with costs `D=2`, `T=3`, and `L=1`, use

> python Diameter.py tree1.newick 2 3 1

To log the result to the end of a log file `Results.csv`, use

> python Diameter.py -l Results tree1.newick 2 3 1 

To calculate the diameters for files `tree1.newick`, `tree2.newick`, and `tree4.newick`, use

>python Diameter.py -r 1 5 tree#.newick 2 3 1

The program will skip over the non-existent `tree3.newick`.

To calculate the diameters for every file in the TreeLifeData folder, logging the result to COG_results.csv, while suppressing text and showing tables, use

>python Diameter.py -r 1 -l COG_results -qd 6000 TreeLifeData/COG####.newick 2 3 1

### Via Interactive Mode:

#### Single Files

To calculate the diameter for a single file in interactive mode, call the following function:
> calculate_diameter_from_file(filename, D, T, L, log=None, debug=False, verbose=True)

The arguments in this function are similar to the command line option `calc`. `filename` is the path to the file that will be reconciled by `DTLReconGraph` and the DTL Diameter calculated for. `D`, `T`, and `L` are the event costs, and `log` is the a string containing path to the csv file used for logging (use `None` for no logging).

`debug`, when set to true, will make Diameter print out the dynamic programming tables it uses to find the solution (but only the ones that are less than 30x30 size).

`verbose`, when set to true, will make Diameter print textual output about the results and their running times.

#### Many Files

To calculate a set of numbered files, use this function:
> repeatedly_calculate_diameter(file_pattern, start, end, d, t, l, log=None, debug=False, verbose=True)

Where `file_pattern` is the pattern used to find the right files (as described in the command line `-i` flag section), `start` is the starting file number, `end` is the exclusive ending file number, and the rest of the parameters are the same as `calculate_diameter_from_file()`.

