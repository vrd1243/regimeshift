#!/bin/python

from rk4 import rk4
import numpy as np
import matplotlib
import sys
from matplotlib import pyplot as plt
from property import get_windowed_property
import ordinal_TSA as ot
from get_data import *

def main():
    
    property = sys.argv[1];
    data_size = 12000;
    window = 6000;
    window_shift = 200;

    [shift_markers, series] = get_data_hybrid1(data_size = data_size, remove_transient = True, normalized = True);
    np.savetxt('series_' + property + '.txt', series); 

    print("Series generated of size ", series.shape);
    
    plt.figure();

    plt.plot(series[:1000], label = 'series');
    plt.legend();
    plt.savefig('real_time_series.png');

    for w in range(500, 8000, 500):

        print(w);
        [x_window, tau] = get_windowed_property(series, w, window_shift, property);

        plt.figure();
        plt.plot(x_window, tau, label='window {}'.format(w));
        plt.title('Window = {}'.format(w));
        for s in shift_markers:
            plt.axvline(x=s, color='g');
        plt.legend();
        plt.savefig('real_time_{}_{}.png'.format(property, w));

main();
