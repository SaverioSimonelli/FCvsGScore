FCvsGScore Created by Saverio Simonelli MIT license 2020

FCvsGScore Python library consists of three main functions that help to merge the information provided from an experiment with the Quad 
(https://github.com/SaverioSimonelli/Quad) output and to graphically represent two variables: y-variable, derived from the experiment, and G-Quadruplexes score.
The three functions are:
Complete Information, Point Graph, Color Graph.

<pre>
Be
t table with quadruplexes information
τ table with differentially expressed transcripts information
r generic row on table t
κ generic information about transcripts on table τ
t' table with quadruplexes and differential expression information

Algorithm: COMPLETE INFORMATION
INPUT: t, τ
OUTPUT: t'
1. set t' := ∅
2. for each row r ∈ t do
3.	get κ in τ relative to r
4. 	replace t' by t' ∪ { κ x r }
5. end for
6. return t'

############

Be
a algorithm
t table with quadruplexes and genes information
x set of quadruplexes scores in a
y set of gene characters in a
G graph

Algorithm: POINT GRAPH
INPUT: a, t
OUTPUT: G
1. get scores x from t about a
2. get characters y from t about a
3. add x and y on G
4. return G

###########

Be
a algorithm
t table with quadruplexes and genes information
x set of quadruplexes scores in a
y set of gene characters in a
C set of colors
G graph

Algorithm: COLOR GRAPH
INPUT: a, t
OUTPUT: G
1. create set C from t
2. get scores x from t about a
3. get characters y from t about a
4. add x and y and C on G
5. return G

</pre>
