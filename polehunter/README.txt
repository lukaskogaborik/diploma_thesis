Short manual of snarkhunter
---------------------------

The latest version of snarkhunter can be obtained from http://caagt.ugent.be/cubic/

Snarkhunter is written for the efficient generation of all nonisomorphic connected simple cubic graphs, but it can be restricted to generating only graphs that have a fixed minimal girth or have certain connectivity requirements.
It can also be used to generate snarks efficiently.

This program has been tested on Linux and Mac OS X.


Installation
------------
First download, extract and configure nauty from http://cs.anu.edu.au/~bdm/nauty/
Snarkhunter requires nauty 2.7 (or a more recent version).

Compile the nauty libraries using: "make nautyW1.a" and "make nautyL1.a"

Then copy the following files to the snarkhunter directory:
nausparse.h nauty.h nautyW1.a nautyL1.a

You can make a snarkhunter binary using "make".
If you want to generate graphs with more than 32 vertices, make a binary using "make 64bit".
This "64-bit binary" can also be used to generate graphs with <= 32 vertices, but is (slightly) slower than the regular binary.


Usage
-----
An overview of all options can also be found by executing "./snarkhunter -h".


Usage: ./snarkhunter <number_of_vertices> <min_girth> [options]

At least the number of vertices and the minimum girth must be specified.
So "./snarkhunter 20 4" generates all cubic with 4,6,...,20 vertices which have girth at least 4.
By default the generated graphs are written to files with names such as: Generated_graphs.20.04 (i.e. cubic graphs with 20 vertices and girth at least 4) in multicode format.
The graph encoding formats are described in formats.txt


Possible options are:

  S: Only output class 2 graphs. This are graphs which are not 3-edge-colourable (i.e. graphs with chromatic index 4).

  C<x>: Only output graphs with connectivity at least x.
        If x >= 4, only graphs with cyclic edge-connectivity at least x are output.

  s: Only output the graphs with 'n' vertices (instead of all graphs with <= n vertices).

  o: The graphs are written to stdout instead of to files.

  g: The graphs will be output in graph6 format instead of multicode format.

  n: No graphs will be output (i.e. the graphs are only counted).

  p: Prime graphs will also be output (to a separate file).

  m <rest> <modulo>: Splits the generation in <modulo> (more or less equally big) parts. Here part <rest> will be executed. 
                     Splitting the generation causes a small overhead, so the sum of the timings for the small parts will be slightly more than the time needed to run the same case without modulo. But this overhead is usually negligible compared to the total execution time.
                     The normal rules for modulo calculation apply. So m 0 2 will give the same result as m 0 4 and m 2 4 combined.


Examples
--------

"./snarkhunter 22 4 S s C4"
Generates all class 2 cubic graphs with 22 vertices and girth at least 4 which are cyclically 4-edge connected.


"./snarkhunter 22 5 o g"
Generates all cubic graphs with at most 22 vertices and outputs them to stdout in graph6 format.


"./snarkhunter 32 4 S s C3 m 10 200"
Splits the generation of all 3-connected class 2 cubic graphs with 32 vertices and girth at least 4 into 200 (more or less equally big) parts and executes part 10.


Don't hesitate to contact me at jan.goedgebeur@ugent.be if you have any further questions or suggestions.


Plugins
-------

If you only wish a fraction of the output, you can save i/o time by writing a plugin. An example is provided in the file edgepent.c, which filters out any graph which contains an edge not lying on any cycle of length at most 5.

The plugin should define a function (or parameter-less macro for the function name) FILTER with a prototype like   int FILTER(GRAPH g, int n).  n is the number of vertices.  g is a 2-dimensional array containing the graph. Vertices are numbered from 0. The neighbours of vertex are g[i][0], g[i][1] and g[i][2], which have type unsigned char.

FILTER is called with each completed cubic graph. If you want to keep this graph, return a non-zero value. If you want to reject it, return 0.

You can optionally provide a function void SUMMARY() which is called just before the final messages are written by snarkhunter. You can use it to write information that was collected by FILTER.

See the makefile for how to compile the plugin into snarkhunter. Note that the double quotation marks have to be protected from the shell. You can use -DPLUGIN='"edgepent.c"' or -DPLUGIN=\"edgepent.c\" in most shells.


Changelog
---------

2025-02-12: Added facility for plugin (kindly added by Brendan McKay).

2015-11-26: Added specialised algorithms for generating cubic graphs and snarks with girth at least 6 or 7.

2011-06-17: First release.
