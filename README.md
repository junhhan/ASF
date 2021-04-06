# Estimation of parameters related to ASF spread in wild boars
### Introduction
The purpose of the script is to explore the posterior distributions of ecological and epidemiological parameters, such as:
<ol><li> transmission rates for infectious live wild boars (WB) or WB carcasses, and</li></ol>
<ol><li> efficiency of roads, rivers, and fences to preventing the spread of ASF,</li></ol>
in the WB population of Gyeonggi and Gangwon provinces, South Korea.

The parameters were estimated with a Bayesian approach. Specifically, we used a Markov Chain Monte Carlo (MCMC) method with a Metropolis-Hastings (MH) algorithm.
<br/><br/>

### How to run
The simulation model was coded in the C programming language. Compilation of .c file requires "mtwister.h" header file, which contains Mersenne Twister pseudo-random number generator based on http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html. Once compiled, one can run the code with a generated executable file.

The executable requires multiple CSV files:
<ol> (1) grid_info.csv, containing cell-wise information,</ol>
<ol> (2) adj_cid.csv, containing the cell ids of close cells,</ol>
<ol> (3) grid_xy.csv, contating the centroid coordinates for each cell,</ol>
<ol> (4) road_cross.csv, containing whether two cells are separated by roads,</ol>
<ol> (5) river_cross.csv, containing whether two cells are separated by rivers,</ol>
<ol> (6) SS_n.csv, containing the number of collected WB carcass samples from a given cell in a given week,</ol>
<ol> (7) SS_k.csv, containing the number of collected ASF (+) WB carcass samples from a given cell in a given week,</ol>
<ol> (8) SS_n_h.csv, containing the number of collected live WB samples from a given cell in a given week,</ol>
<ol> (9) SS_k_h.csv, containing the number of collected ASF (+) live WB samples from a given cell in a given week.</ol>
Among them, the files in (6) ~ (9) cannot be freely distributed due to the lack of authorisation.
<br/><br/>

### Notes
The C code was designed to use parallel computation with 2 cores to acheive the computational efficiency. 
Also, the model was coded to iterate the simulation for 70000 times, which requires approximately __*72 hours*__.
