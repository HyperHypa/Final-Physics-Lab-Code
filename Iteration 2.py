import numpy as np
import matplotlib.pyplot as plt

# https://medium.com/modern-physics/simple-pendulum-odesolver-using-python-dcb30c267eee

class ODESolver(object):
    """Second-order ODE Solver.
    Parameters
    ------------
    omega_0 : float
            initial angular velocity
    theta_0 : float
            initial angular displacement
    eta : float
        time step size
    n_iter : int
           number of steps
        
    Attributes
    -----------
    time_ : 1d-array
        Stores time values for each time step.
    omega_ : 1d-array
        Stores angular velocity values for each time step.
    theta_ : 1d-arra
       Stores angular displacement values for each time step.
        
    Methods
    -----------
    euler(alpha): Implements the Euler algorithm for the acceleration function alpha.
    
    midpoint(alpha): Implements the Midpoint algorithm for the acceleration function alpha.
    
    verlet(alpha): Implements the Verlet algorithm for the acceleration function alpha.
    """
    def __init__(self, omega_0, theta_0, phi_0, eta, n_iter, theta_mass, phi_mass, spring_distance_theta, spring_distance_phi, spring_constant, spring_length):
        self.omega_0 = omega_0
        self.theta_0 = theta_0
        self.phi_0 = phi_0
        self.eta = eta
        self.n_iter = n_iter
        self.theta_mass = theta_mass
        self.phi_mass = phi_mass
        self.spring_distance_theta = spring_distance_theta
        self.spring_distance_phi = spring_distance_phi
        self.spring_constant = spring_constant
        self.spring_length = spring_length
        
    
    
    def verlet(self,alpha):
        """Implement Verlet Method.
        
        Parameters
        ----------
        alpha : acceleration function
        Returns
        -------
        self : object
        """
        self.time_ = np.zeros(self.n_iter)
        self.theta_ = np.zeros(self.n_iter)
        self.phi_ = np.zeros(self.n_iter)
        self.theta_[0] = self.theta_0*np.pi/180.0
        self.phi_[0] = self.phi_0*np.pi/180.0
        self.time_[1]= self.eta
        self.theta_[1] = self.theta_[0]+self.omega_0*self.eta +0.5* (self.eta**2)*alpha(self.theta_[0], self.phi_[0], self.spring_length, self.spring_constant, self.theta_mass)
        self.phi_[1] = self.phi_[0]+self.omega_0*self.eta +0.5* (self.eta**2)*alpha(self.phi_[0], self.theta_[0], self.spring_length, self.spring_constant, self.phi_mass)
        
        for i in range(self.n_iter-2):
            self.time_[i+2] = self.time_[i+1] + self.eta
            self.theta_[i+2] = 2.0*self.theta_[i+1] -self.theta_[i] + (self.eta**2)*alpha(self.theta_[i+1], self.phi_[i+1], self.spring_length, self.spring_constant, self.theta_mass)
            self.phi_[i+2] = 2.0*self.phi_[i+1] -self.phi_[i] + (self.eta**2)*alpha(self.phi_[i+1], self.theta_[i+1], self.spring_length, self.spring_constant, self.phi_mass)
        return self

def alpha(x, y, sl, sc, m):
    g = 9.80665 # m/s^2
    l = 10 # m
    dis = ((np.cos(x) - np.cos(y))**2 + (np.sin(x) - np.sin(y) + sl)**2)**0.5
    return -g / l * np.sin(x) - (dis - sl)*sc/m * np.sin((np.sin(x) - np.sin(y) + sl)/dis)

pendulum = ODESolver(omega_0 = 0, theta_0 = 10, phi_0 = 0, eta=0.1, n_iter=100, theta_mass = 1, phi_mass = 1, spring_distance_theta = 5, spring_distance_phi = 5, spring_constant = 0.5, spring_length = 2).verlet(alpha)
time = pendulum.time_
theta = pendulum.theta_
phi = pendulum.phi_
plt.plot(time,theta*180/np.pi,lw=3,color='blue', label = 'theta')
plt.plot(time,phi*180/np.pi,lw=3,color='red', label = 'phi')
plt.legend()
plt.xlabel('time(s)',size=13)
plt.ylabel('angle (deg)',size=13)
plt.title('Verlet Method',size=13)
plt.show()