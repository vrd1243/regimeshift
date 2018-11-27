#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

int main(int argc, char **argv) {
    
    int i = 0;
     
    int tau_avg = 5;
    int n_neighbors = 3; 
    int length = atoi(argv[3]);
    int tau_max = 100;
    int tau_start = 1;

    double series[length];
    double menger_sum[tau_max];
    double menger_variance[tau_max];

    ifstream f_in(argv[1]);
    ofstream f_out(argv[2]);
    
    while(1) {
	f_in >> series[i++];
	if (f_in.eof()) break;
    }
    
    double m = 0;
    float total_sum = 0;
    float total_variance = 0;

    for (int tau = tau_start; tau < tau_max; tau++) {
      
      int menger_length = length - 2*tau_avg - 2*n_neighbors - tau;
      double menger[menger_length];	  

      for (int t = 0; t < menger_length; t++) {
          menger[t] = 0;	
      }

      for (int t = tau_avg + n_neighbors; t < length - tau - tau_avg - n_neighbors; t++) {

          m = 0; 

	  for (int k = -n_neighbors; k < n_neighbors + 1; k++) {

	      double x_first = series[t-tau_avg+k], y_first = series[t+tau-tau_avg+k];
	      double x_second = series[t+k], y_second = series[t+tau+k];
	      double x_third = series[t+tau_avg+k], y_third = series[t+tau+tau_avg+k];
	      
	      double a = sqrt((x_first - x_second)*(x_first - x_second) + 
			     (y_first - y_second)*(y_first - y_second));

	      double b = sqrt((x_third - x_second)*(x_third - x_second) + 
			     (y_third - y_second)*(y_third - y_second));

	      double c = sqrt((x_first - x_third)*(x_first - x_third) + 
			     (y_first - y_third)*(y_first - y_third));

	      double s = (a + b + c)/2;

	      //if (isnan(sqrt(s*(s-a)*(s-b)*(s-c))))
	      //	cout << a << " " << b << " " << c << " " << s << " " << s*(s-a)*(s-b)*(s-c) << endl;
	      
	      if (((a*b*c) != 0) && (s*(s-a)*(s-b)*(s-c) >= 0))
		  m += 4*sqrt(s*(s-a)*(s-b)*(s-c)) / (a*b*c);
	     
	      //cout <<  4*sqrt(s*(s-a)*(s-b)*(s-c)) / (a*b*c) << endl;
	  }

	  menger[t - tau_avg - n_neighbors] = m / (2*n_neighbors + 1);
      }
      
      menger_sum[tau] = 0;
     
      for (int t = 0; t < menger_length; t++) {
          menger_sum[tau] += menger[t];	
      }

      menger_variance[tau] = 0;
      
      double menger_mean = menger_sum[tau] / (menger_length);

      for (int t = 0; t <  menger_length; t++) {
          menger_variance[tau] += (menger[t] - menger_mean)*(menger[t] - menger_mean);	
      }

      menger_variance[tau] /= menger_length;
      
      total_sum += menger_sum[tau];
      total_variance += menger_variance[tau];

    }

    for(int tau = tau_start; tau < tau_max; tau++) {
        f_out << menger_sum[tau] / menger_sum[tau_start] << " " << menger_variance[tau] / menger_variance[tau_start] << endl;  
    }
       
}
