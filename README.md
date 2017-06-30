# dtldiameter

## __**How to use DTLReconGraph.py**__
One could run this file in interactive mode and test any functions individually; however, the most computationally important function in this file is DP, especially since its results are used in Diameter.py. The most helpful function to run, however, is reconcile - this is what Diameter.py actually uses to compute results (including the DTL reconciliation graph itself and the number of MPRs for a given data set), and it extracts all of the important info from DP and allows the user to pass data files for the program to run on. Consequently, users should focus on understanding the reconcile function and all functions on which it depends in order to get the most out of this program. Of course, all of the other functions that don't directly have anything to do with reconcile would be useful to understand also. All inputs, return values, and their formats (along with the rest of the code) are thoroughly explained in the file in the form of comments and docstrings.

## __**How to use Diameter.py**__

### From the command line:

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

### Via interactive mode:
