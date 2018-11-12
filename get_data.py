#!/bin/python

from rk4 import rk4
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from datasets.lorenz import get_lorenz, set_params
import ordinal_TSA as ot
from datasets.pendulum import derivative, changeA
from datasets.rossler import get_rossler

def normalize(series, normalized = False):
    
    if normalized:
        return ((series - np.mean(series)) / np.std(series));
    else: 
        return series;

def get_data_lorenz(data_size, normalized = False):
    [t, data] = get_lorenz(n=data_size, tmax=500);
    plt.figure();
    plt.plot(data[:,0], data[:,1]);
    plt.savefig('lorenz.png');
    np.savetxt('lorenz.txt', data[:,0]);
    return [[], data[:,0]];

def get_data_hybrid1(data_size, remove_transient = False, normalized = False):    
    
    shift_markers = [];
    
    if remove_transient:
        start = int(data_size / 6);
    else:
        start = 1;

    [t, data] = get_lorenz(n=data_size, tmax=500);
    series = normalize(data[start:,0], normalized);
    end_point = data[-1,:2];
    plt.figure();
    plt.plot(data[start:,0], data[start:,1], '.')
    plt.savefig('hybrid1_lorenz1.png');
    np.savetxt('hybrid1_lorenz1.txt', data[start:,0]);

    shift_markers.append(len(series));

    A_list = np.arange(0.86,0.87,0.02)
    for A in A_list:
        print("A=", A);
        changeA(A);
        [data1,t] = rk4(derivative, end_point, 0, 0.005, data_size); 
        plt.figure();
        plt.plot(data1[start:,0], data1[start:,1], '.')
        plt.savefig('hybrid1_pendulum.png');
        np.savetxt('hybrid1_pendulum.txt', data[:,0]);
        series1 = normalize(data1[start:,0], normalized);
        series = np.concatenate((series, series1), axis=0);
        shift_markers.append(len(series));
        end_point = data1[-1,:];

    [t, data] = get_lorenz(n=data_size, tmax=500, init = [end_point[0], end_point[1], 0]);
    series1 = normalize(data[start:,0], normalized);
    series = np.concatenate((series, series1), axis=0);

    plt.figure();
    plt.plot(data[start:,0], data[start:,1], '.')
    plt.savefig('hybrid1_lorenz2.png');
    np.savetxt('hybrid1_lorenz2.txt', data[start:,0]);
     
    end_point = data1[-1,:];
    shift_markers.append(len(series));

    [t, data] = get_rossler(n=data_size, tmax=500, 
                start = [end_point[0], end_point[1], 0]);
    
    plt.figure();
    plt.plot(data[start:,0], data[start:,1], '.')
    plt.savefig('hybrid1_rossler_2d.png');
    np.savetxt('hybrid1_rossler.txt', data[start:,0]);
    
    series1 = normalize(data[start:,0], normalized);
    series = np.concatenate((series, series1), axis=0);

    return [shift_markers, series];

def get_data_hybrid2(data_size, remove_transient = False, normalized = False):
    
    shift_markers = [];
    
    end_point = [1,3.14];

    if remove_transient:
        start = int(data_size / 6);
    else:
        start = 1;

    changeA(0.47);
    [data,t] = rk4(derivative, end_point, 0, 0.005, data_size); 
    series = normalize(data[start:,0], normalized);
    shift_markers.append(len(series));
    end_point = data[-1,:];
    plt.figure();
    plt.plot(data[start:,0], data[start:,1], '.');
    plt.savefig('pendulum_A_{}.png'.format(0.47));

    A_list = [0.6, 0.86, 0.89, 0.92]
    for A in A_list:
        print("A=", A);
        changeA(A);
        [data,t] = rk4(derivative, end_point, t[-1], 0.005, data_size); 
        series1 = normalize(data[start:,0], normalized);
        series = np.concatenate((series, series1), axis=0);
        shift_markers.append(len(series));
        end_point = data[-1,:];
        plt.figure();
        plt.plot(data[start:,0], data[start:,1], '.');
        plt.savefig('pendulum_A_{}.png'.format(A));

    #changeA(0.86); 
    #[data,t] = rk4(derivative, [3.14, 50], 0, 0.005, data_size); 
    #series = np.concatenate((series, data[:,0]), axis=0);
    #series = data[:,0]

    return [shift_markers, series];

def get_data_pendulum(data_size, remove_transient = False, normalized = False):

    if remove_transient:
        start = int(data_size / 6);
    else:
        start = 1;

    shift_markers = [];
    
    end_point = [3.14,50]; 

    changeA(0.91);
    [data,t] = rk4(derivative, end_point, 0, 0.005, data_size); 
    series = normalize(data[start:,0], normalized);
    end_point = data[-1,:];
    shift_markers.append(len(series));

    plt.figure();
    plt.plot(data[start:,0], data[start:,1], '.');
    plt.savefig('pendulum_A_{}.png'.format(0.91));
    
    A_list = [.6, 0.886];

    for A in A_list:
        print("A=", A);
        changeA(A);
        [data,t] = rk4(derivative, end_point, t[-1], 0.005, data_size); 

        plt.figure();
        plt.plot(data[start:,0], data[start:,1], '.');
        plt.savefig('pendulum_A_{}.png'.format(A));

        series1 = normalize(data[start:,0], normalized);
        series = np.concatenate((series, series1), axis=0);
        shift_markers.append(len(series));
        end_point = data[-1,:];
    
    np.savetxt('pendulum_hybrid.txt', series);

    return [shift_markers, series]

def get_data_sinusoidal(data_size, remove_transient = False, normalized = False):

    
    if remove_transient:
        start = int(data_size / 6);
    else:
        start = 1;

    shift_markers = [];
    
    t = np.arange(0,1000,0.1);
    
    first = np.sin(t); 
    second = np.sin(t) + np.sin(2*t);
    third = np.sin(t) + np.sin(2*t) + np.sin(3*t);
    
    np.savetxt('sine_first.txt', first);
    np.savetxt('sine_second.txt', second);
    np.savetxt('sine_third.txt', third);

    series = first;
    shift_markers.append(len(series));
    series = np.concatenate((series, second), axis = 0);
    shift_markers.append(len(series));
    series = np.concatenate((series, third), axis = 0);
    shift_markers.append(len(series));
    
    return [shift_markers, series]
