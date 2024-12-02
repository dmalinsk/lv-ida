# lv-ida
For the implementation of LV-IDA in R
@author: Daniel Malinsky d.malinsky@columbia.edu

 "LV-IDA for estimating causal effects with ancestral graph Markov Models."
 Copyright (C) 2016  Daniel Malinsky

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


This repository has the code for LV-IDA (Latent Variable IDA). The algorithm is described in:

Daniel Malinsky and Peter Spirtes. “Estimating bounds on causal effects in high-dimensional and possibly confounded systems,” 
International Journal of Approximate Reasoning 88: 371-384, 2017.

Daniel Malinsky and Peter Spirtes. "Estimating causal effects with ancestral graph Markov models." 
Journal of Machine Learning Research: Workshop and Conference Proceedings (PGM) 52, 2016. 
(Proceedings of the Eighth International Conference on Probabilistic Graphical Models)

It is based on IDA by Maathuis, Kalisch, and Buhlmann. As far as the end-user is concerned, it works just like IDA as implemented in the R package pcalg. 
See the paper by Kalisch et al. "Causal inference using graphical models with the R package pcalg." Journal of Statistical Software, 47(11):1-26, 2012. 
If you follow the instructions for using IDA from that paper, you can use LV-IDA.

The important differences:

1) LV-IDA takes a PAG (Partial Ancestral Graph) as input, e.g., the output of the FCI algorithm. 
(IDA takes a CPDAG, e.g., the output of the PC or GES algorithms.)

2) LV-IDA takes an adjacency matrix instead of a graph as the input representation of a PAG. So, if you run:
> fci.est <- fci( ...arguments... )

then do:

> lv.ida(x,y,cov(data),fci.est@amat,method="local")

to estimate the multiset of total causal effects of x on y. fci.est@amat gets the right adjacency matrix. 
x and y should be numbers corresponding to vertices in the graphical model. 
(Note: the columns and rows of the adjacency matrix should be numbers, in order. 
This version does not currently support adjacency matrices with character variable names as row/column names.)

If you're not using FCI or RFCI or some other PAG search method in pcalg, then you might need to turn your graph 
search result into the right kind of adjacency matrix. 1/2/3 at (i,j) means circle/arrowhead/tail at j from i. 0 means no adjacency.
For example:

graph: 1<-2<-o3--4

matrix:
  1 2 3 4 \n
1 0 3 0 0 \n
2 2 0 1 0 \n
3 0 2 0 3 \n
4 0 0 3 0 \n

No cycles are allowed! Sometimes search algorithms will return cyclic graphs even when that violates the assumption of an acyclic generating structure. 
The function "iscyclic" in iscyclic.R can be used to check that you have not got a cyclic graph.

Try running the code in example.R to get a feel for how LV-IDA works. Write d.malinsky@columbia.edu with questions/problems.

NOTE (Dec 1 2024): an important bug was fixed in late 2024. Among other things, there was a problem in how the "local" version of LV-IDA was implemented, the local set described in the published paper was not actually sufficient to check all the relevant graphical conditions for the backdoor criterion. The update makes this set somewhat larger, which could affect computation time. An erratum describing this change more formally is in development and this page will be updated when that is posted.
