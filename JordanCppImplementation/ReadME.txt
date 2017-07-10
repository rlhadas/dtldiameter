"tree.cpp" conatins pretty much everything, from the representations
of trees to the dp algorithm and reconciliation graphs/trees.
It also has a parser for .newick files. 

Please note that the overall choices for representing R-Graphs as
a 2D vector of indexed mapping nodes was made in the interest of speed.
It is quite a bit faster than a map.


"main.cpp" runs all of the statistics over the TreeLifeData set. It
outputs a .csv format that you can read in excel. You'll have to 
hardcode in any file changes.

"draw.cpp" includes tools for drawing trees, reconciliation graphs, and
reconciliation trees. It isn't polished. But, it may be helpful
if you need to visualize something. You might need a specific c++
library to run it (maybe).