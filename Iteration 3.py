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
        initial angular displacement of the right pendulum
    phi_0 : float
        initial angular displacement of the left pendulum
    eta : float
        time step size
    n_iter : int
        number of steps
    theta_mass : float
        mass of the right pendulum's bob
    phi_mass : float
        mass of the left pendulum's bob
    system_mass : float
        mass of the the system excluding the bobs

        
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
    verlet(alpha): Implements the Verlet algorithm for the acceleration function alpha.
    """
    def __init__(self, omega_0, theta_0, phi_0, eta, n_iter, theta_mass, phi_mass):
        self.omega_0 = omega_0
        self.theta_0 = theta_0
        self.phi_0 = phi_0
        self.eta = eta
        self.n_iter = n_iter
        self.theta_mass = theta_mass
        self.phi_mass = phi_mass
        
    
    
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
        self.theta_[1] = self.theta_[0]+self.omega_0*self.eta +0.5* (self.eta**2)*alpha(self.theta_[0], self.phi_[0], self.theta_mass, self.phi_mass)
        self.phi_[1] = self.phi_[0]+self.omega_0*self.eta +0.5* (self.eta**2)*alpha(self.phi_[0], self.theta_[0], self.phi_mass, self.theta_mass)
        
        for i in range(self.n_iter-2):
            self.time_[i+2] = self.time_[i+1] + self.eta
            self.theta_[i+2] = 2.0*self.theta_[i+1] -self.theta_[i] + (self.eta**2)*alpha(self.theta_[i+1], self.phi_[i+1], self.theta_mass, self.phi_mass)
            self.phi_[i+2] = 2.0*self.phi_[i+1] -self.phi_[i] + (self.eta**2)*alpha(self.phi_[i+1], self.theta_[i+1], self.phi_mass, self.theta_mass)
        return self

def alpha(x, y, mx, my):
    g = 9.80665 # m/s^2
    l = 10 # m
    asq = (-g/l*np.sin(x) * np.sin(x) + -g/l*np.sin(y) * np.cos(y) * mx / my)**2 + -g/l*np.sin(x) * np.cos(x)
    if asq > 0:
        a = asq**0.5
    else:
        a = -(-asq)**0.5
    print(x)
    return a

pendulum = ODESolver(omega_0 = 0, theta_0 = 10, phi_0 = 0, eta=0.1, n_iter=100, theta_mass = 1, phi_mass = 1).verlet(alpha)
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