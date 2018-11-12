#!/bin/python

import numpy as np
import ordinal_TSA
import matplotlib
from matplotlib import pyplot as plt
plt.switch_backend('agg');
from pytisean import tiseanio;
from curvature import get_menger_centered_avg 
import numpy as np
import time
from subprocess import call

# All those million thresholds!
stable_tau_threshold = 0.1;
minimum_threshold = 0.05;
sufficient_window_threshold = 500;
max_window_size = 1000;
window_shift = 20;

def get_windowed_property(series, window, window_shift, property): 
    
    tau = [];
    tau_max = 100;

    if property == "mi":
        get_best_property = get_mi_tau;
    elif property == "menger":
        get_best_property = get_menger_tau;
    
    x_window = np.arange(0,len(series) - window, window_shift)
    for i in x_window:
        full = series[i:i+window];
        tau.append(get_best_property(full, tau_max));

    return [x_window + window , tau];

def find_minimum(property):
  
    x_prev = property[0]; 
    for idx, x_curr in enumerate(property[1:]):
        if x_curr < x_prev:
            break;
        
        x_prev = x_curr; 
     
    x_first = property[idx];
    x_min = x_first;
    first_idx = idx;
    min_idx = first_idx;

    for idx, x_curr in enumerate(property[first_idx + 1:]):
        
        if x_curr > x_min and x_first == x_min:
            print(x_first, first_idx, idx, x_curr, x_min);
            print(property);

        if x_curr > x_min and (abs(x_min - x_curr) / abs(x_first - x_min) >= minimum_threshold):
            break;

        if x_curr <= x_min:
            x_min = x_curr;
            min_idx = idx + 1;
    
    min_idx += first_idx + 1;
    
    return min_idx + 1;

def get_mi_tau(data, tau_max = 100):
    mutual, err = tiseanio('mutual', '-D', tau_max, data=data, silent=True);
    print(mutual.shape, err)
    return (find_minimum(mutual[:,1]));
    
def get_menger_tau(data, tau_max):
    menger = [];
    np.savetxt('menger_in.txt', data);
    call(["./menger", "menger_in.txt", "menger_out.txt", str(data.shape[0])])
    menger_data = np.loadtxt("menger_out.txt");

    return (find_minimum(menger_data[:,0]));

def get_window(data, start):

    tau_list = [];

    tau_prev = -np.inf;
    tau_start = -np.inf;
    counter = 0;
    tau_max = 100;

    if len(data) - start < max_window_size:
        max_window = len(data)- start;
    else:
        max_window = max_window_size;
    
    window = 0;
    for window in np.arange(tau_max, max_window, window_shift):
        tau = get_mi_tau(data[start:start+window], tau_max); #(find_minimum(mutual[:,1]));
        tau_list.append(tau);
        if tau_max < tau:
            tau_max = tau;
        
        if np.abs(tau - tau_start) < stable_tau_threshold * tau_max:
            counter += window_shift;

        else:
            if counter >= sufficient_window_threshold:
                window -= window_shift;
                break;
            counter = 0;
            tau_start = tau;
    
        tau_prev = tau;
    
    if window != 0: 
        plt.figure();
        plt.plot(mutual[:,1]);
        plt.savefig('tau_plot_%d_%d.png' % (start, start + window));

    #plt.figure();
    #plt.plot(tau_list);
    #plt.savefig('tau_plot_%d_%d.png' % (start, start + window));
    #time.sleep(2);
    return (tau_start, window);

def embed_window(series, tau, start, window):
    embed, err = tiseanio('delay', '-m', 2, '-d', tau, data=series, silent=True);
    return embed;
