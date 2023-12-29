import numpy as np
import matplotlib.pyplot as plt


# https://medium.com/modern-physics/simple-pendulum-odesolver-using-python-dcb30c267eee

def equation(x, y, f, k):
    return f ** 0.5 * 2 * np.pi + k * np.sin(y - x)


class Pendulum(object):
    """
    Parameters
    ------------
    theta_0 : float
        initial angular displacement of the left pendulum
    phi_0 : float
        initial angular displacement of the right pendulum
    eta : float
        time step size
    n_iter : int
        number of steps
    frequency:
        frequency of both pendulums
    coupling_constant:
        the coupling strength of the pendulum

        
    Attributes
    -----------
    times : 1d-array
        Stores time values for each time step.
    thetas : 1d-array
       Stores angular displacement values for each time step of the left pendulum.
    phis : 1d-array
       Stores angular displacement values for each time step of the right pendulum.
        
    Methods
    -----------
    kuramoto_model:
        Runs the Kuramoto model.
    kuramoto:
        The equation of the Kuramoto model.
    """

    def __init__(self, theta_0, phi_0, eta, n_iter, frequency, coupling_constant):
        self.theta_0 = theta_0
        self.phi_0 = phi_0
        self.eta = eta
        self.n_iter = n_iter
        self.frequency = frequency
        self.coupling_constant = coupling_constant
        self.times = np.zeros(self.n_iter)
        self.thetas = np.zeros(self.n_iter)
        self.phis = np.zeros(self.n_iter)

    def kuramoto_model(self):
        self.thetas[0] = self.theta_0 * np.pi / 180.0
        self.phis[0] = self.phi_0 * np.pi / 180.0

        for i in range(self.n_iter - 1):
            self.times[i + 1] = self.times[i] + self.eta
            self.thetas[i + 1] = self.thetas[i] + self.eta * equation(self.thetas[i], self.phis[i], self.frequency,
                                                                      self.coupling_constant)
            # print(self.equation(self.thetas[i], self.phis[i], self.theta_length, self.coupling_constant))
            self.phis[i + 1] = self.phis[i] + self.eta * equation(self.phis[i], self.thetas[i], self.frequency,
                                                                  self.coupling_constant)
            # print(self.equation(self.phis[i], self.thetas[i], self.phi_length, self.coupling_constant))
        return self


system = Pendulum(theta_0=-90, phi_0=89, eta=0.001, n_iter=700, frequency=184*60/2, coupling_constant=6.5).kuramoto_model()
time = system.times * 180 / np.pi
left_pendulum = np.sin(system.thetas) * 0.127
right_pendulum = np.sin(system.phis) * 0.127

plt.figure(figsize=(14, 6.4))
plt.plot(time, left_pendulum, lw=3, color='blue', label='Left pendulum')
plt.plot(time, right_pendulum, lw=3, color='red', label='Right pendulum')
plt.legend()
plt.xlabel('Time (s)', size=13)
plt.ylabel('X-Position (m)', size=13)
plt.title('X-Position Against Time', size=13)
plt.show()
