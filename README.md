# Chapman-Kolmogorov-test
This is a Python script to do Chapman-Kolmogorov test in Markov State Model building of molecular dynamics using MSMBuilder.

It uses the transition probability between different states (tProb.mtx), population of each Markov state (Populations.dat) and the assignment of snapshots to different states (Assignments.Fixed.h5) to test whether the model fits the data. It plots the remaining probability after certain number of Markov steps based on the model and and also that based on the data with error estimation. The users need to determine the proper block size used in block averaging for error estimation in highlighted position of the code. The script takes four arguments:
1, number of top populated states to do CK test
2, number of markov steps to do CK test
3, number of primary steps or conformational snapshots to reach Markovianity
4, total number of states in the model

Then the code can be executed by: python ck_test.py #1 #2 #3 #4

Author: Yuqing Zheng (syqzheng@gmail.com)
Jun. 24th, 2016
