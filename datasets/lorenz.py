import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Lorenz paramters and initial conditions
sigma, beta, rho = 10, 8/3, 28
u0, v0, w0 = 0, 1, 1.05

def lorenz(X, t, sigma, beta, rho):
    """The Lorenz equations."""
    u, v, w = X
    up = -sigma*(u - v)
    vp = rho*u - v - u*w
    wp = -beta*w + u*v
    return up, vp, wp

def set_params(s, b, r):
    global sigma, beta, rho
    sigma = s;
    beta = b;
    rho = r;

def get_lorenz(tmax=100, n=10000, init=[0,1,1.05]):
# Integrate the Lorenz equations on the time grid t

  u0 = init[0];
  v0 = init[1];
  w0 = init[2];

  t = np.linspace(0, tmax, n+1)
  f = odeint(lorenz, (u0, v0, w0), t, args=(sigma, beta, rho))
  x, y, z = f.T

  array = np.concatenate((np.reshape(x,(-1,1)), np.reshape(y,(-1,1)),np.reshape(z,(-1,1))), axis=1);
  np.savetxt('lorenz.txt', array);
  return [np.reshape(t,-1,1), array];
