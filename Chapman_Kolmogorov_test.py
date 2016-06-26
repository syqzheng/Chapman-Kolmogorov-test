# This is a script to run Chapman-Kolmogorov test of Markov State Models of molecular dynamics.
# This script reads the output files from MSMBuilder (http://msmbuilder.org/). It uses the
# transition probability between different states, population of each Markov state and the assignment 
# of snapshots to different states to test whether the model fits the data. It plots the remaining
# probability after certain number of Markov steps based on the model and and also that based on the
# data with error estimation.
#
# Author: Yuqing Zheng (syqzheng@gmail.com)
# Jun. 24th, 2016
#

import numpy as np
import matplotlib.pyplot as plt
import h5py
import tables
import sys

def remaining_probability_from_model(num_states, num_steps, total_states, tprob, pop_sort):
    """
    Calculate the remaining probability for top populated states after certain
    number of Markov steps from the Markov State Model.
    """
    prob_model = []
    for i in range(1,num_states+1):
        prob = np.zeros(total_states)
        prob[pop_sort[-i]] = 1.0
        prob_state = []
        for j in range(1, num_steps+1):
            prob = np.dot(prob, tprob)
            prob_state.append(prob[pop_sort[-i]])
        prob_model.append(prob_state)
    return prob_model


def remaining_probability_from_data(num_states, num_steps, markov_step, pop_sort, data):
    """
    Calculate the remaining probability for top populated states after certain
    number of Markov steps from trajectories.
    """
    length = len(data[0])
    state_flag = []
    prob_data = []
    for i in range(1,num_states+1):
        print 'Calculating state', i
        lag_flag = []
        prob_state = []
        for k in range(1,num_steps+1): # to get n x Markov
            markov_step_new = markov_step*k
            state = pop_sort[-i]
            flag = []
            for traj in data:
                for j in range(length - markov_step_new):
                    if traj[j+markov_step_new] != -1 and traj[j] == state:
                        if traj[j+markov_step_new] == state:
                            flag.append(1)
                        else:
                            flag.append(0)
            lag_flag.append(flag)
            prob_state.append(float(sum(flag))/len(flag))
        prob_data.append(prob_state)
        state_flag.append(lag_flag)
    return state_flag, prob_data


def block_average(x, block_size):
    """
    Calculate uncertainty of remaining probability using block averaging.
    """
    block_means = []
    for i in range(0, len(x)/block_size):
        block_means.append(np.mean(x[block_size*i:block_size*(i+1)]))         
    sigma = np.std(block_means)/np.sqrt(len(x)/block_size-1)
    return sigma

def get_standard_error(num_states, num_steps, state_flag, block_sizes):
    """
    Calculate uncertainty of remaining probability for each state at each Markov step.
    """
    standard_error = []
    for i in range(0, num_states):
        error = []
        for j in range(0, num_steps):
            error.append(block_average(state_flag[i][j], block_sizes[i]))
        standard_error.append(error)
    return standard_error


def plot_ck_test(num_states, num_steps, prob_data, standard_error, prob_model):
    """
    Plot the CK test results for each state and save as 'CK_test.png'.
    """
    lag = range(1, num_steps+1)
    plt.figure(figsize=(8,2.5*num_states))
    fs = 10
    
    position = []
    for i in range(1, num_states+1):
        position.append(int(str(num_states)+'1'+str(i)))
    cl1, cl2 = 'blue', 'magenta'
        
    for i in range(0, num_states):
        plt.subplot(position[i])
        plt.plot(lag, prob_data[i], c=cl1, label='From Data', lw=2.5, marker='o', markeredgecolor=cl1, markersize=8)
        plt.errorbar(lag, prob_data[i], c=cl1,  yerr=standard_error[i], ecolor='blue', lw=2.0)
        plt.plot(lag, prob_model[i], c=cl2, label='From Model', lw=2.5, marker='o', markeredgecolor=cl2, markersize=8)
        plt.title('State '+str(i+1), fontsize = fs)
        plt.xlabel('Number of Markov steps', fontsize = fs)
        plt.ylabel('Probability remaining', fontsize = fs)
        plt.tick_params(axis='both', labelsize=fs)
        plt.xlim(0.8, )
        plt.ylim(0.0, 1.0)
    plt.tight_layout()
    plt.savefig('CK_test.png', dpi=600)
    plt.show()

def entry_point():
    num_states = int(sys.argv[1])   # number of top populated states to do CK test
    num_steps = int(sys.argv[2])    # number of markov steps to do CK test
    markov_step = int(sys.argv[3])  # markov_time / time_interval_between_snapshots
    total_states = int(sys.argv[4]) # total number of states
    
    assign = h5py.File('Assignments.Fixed.h5', 'r')
    data = np.array(assign.get("arr_0"))                     # state assignment of each snapshot
    
    pop = np.loadtxt('Populations.dat')                      # Population of each state
    pop_sort = sorted(range(len(pop)), key=lambda k: pop[k]) # state index ranked by population
    
    tprob = np.zeros((total_states, total_states))
    with open('tProb.mtx', 'rb') as f:
        for l in f.readlines()[3:]:
            tprob[int(l.split(' ')[0])-1][int(l.split(' ')[1])-1] = float(l.split(' ')[2])
    # transition probability matrix

    #############################################################
    block_sizes = [1000]*num_states     # Adjust here to give proper block size for each state as a list.
    print 'Using a block size of 1000 for each state, which may not be optimal.'
    print 'Please modify line 124 to use proper block sizes to get proper error estimate.'
    #############################################################

    prob_model = remaining_probability_from_model(num_states, num_steps, total_states, tprob, pop_sort)
    state_flag, prob_data = remaining_probability_from_data(num_states, num_steps, markov_step, pop_sort, data)
    standard_error = get_standard_error(num_states, num_steps, state_flag, block_sizes)
    plot_ck_test(num_states, num_steps, prob_data, standard_error, prob_model)

if __name__ == "__main__":
    entry_point()
