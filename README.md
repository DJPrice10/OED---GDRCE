# OED-GDRCE
Optimal experimental design code for group dose-response challenge experiments (C. jejuni in chickens)

Code folder contains all relevant code to run the different scenarios considered in the manuscript, "Designing group dose-response studies in the presence of transmission".

Note that when evaluating the EKLD, a mex function, "ABCDE\_chunk", is called upon. This file has been compiled for a mac (file ending .mexmaci64), and may need to be compiled if running an alternative OS. The details of this method for evaluating the EKLD can be found in: Price, Bean, Ross, Tuke, _On the efficient determination of optimal Bayesian experimental designs using ABC: A case study in optimal observation of epidemics_, Journal of Statistical Planning & Inference, 172, 1--15, 2016.

Scripts to run scenarios are of the format: INSH\_{prior distribution}\_{target parameters}\_{utility}.mat

Subordinate functions: 
1. ABCDE\_chunk.c: C code to evaluate the EKLD using the ABCdE approach (see above for reference).
2. ABCDE\_chunk.maxmaci64: mex compiled version of above.
3. getdataSEKI.m: Function to call "SEKI\_simulate.m" to simulate the SEKI process via the Gillespie algorithm, and return data for the given observation time(s). Input given _transmission\_rate_, _exposure\_rate_, group size _N_, initial number in first exposed class (i.e., infected from dose) _e1_, maximum time to simulate process up to _Tmax_,  and number of exposed classes _nclasses_.
4. SEKI\_simulate.m: Function containing Gillespie algorithm for the SEKI process. Inputs as per "getdataSEKI.m".
5. infected\_pop\_at\_time.m: Function to return state of process at given observation time(s) _x_, given _state_ and _time_ of transitions from Gillespie algorithm.
6. howmanyunique.m: Function to count the number of unique rows in a matrix _A_.
7. discrete\_trunc\_sample.m: Function to sample _m_ points around current value _T_, across feasible range _minmaxT_, according to a multivariate normal distribution with standard deviation _sd\_T_. Sampled points are rounded to grid according to _stepsize_.



Please also note that the code is correct for the examples we have considered, however, we cannot guarantee that there will not be an issue if aspects of the code are changed to consider alternative input (i.e., other prior distributions). Care must be taken to ensure all aspects of the code are correct for alternative examples.
