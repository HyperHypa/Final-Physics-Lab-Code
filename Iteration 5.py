import numpy as np
import matplotlib.pyplot as plt


# https://medium.com/modern-physics/simple-pendulum-odesolver-using-python-dcb30c267eee

def equation(x, y, length, k):
    g = 9.80665  # m/s^2
    return (2 * g / length) ** 0.5 + k * np.sin(y - x)


class Pendulum(object):
    """
    Parameters
    ------------
    theta_0 : float
        initial angular displacement of the right pendulum
    phi_0 : float
        initial angular displacement of the left pendulum
    eta : float
        time step size
    n_iter : int
        number of steps
    theta_length:
        length of the right pendulum
    phi_length:
        length of the left pendulum
    coupling_constant:
        the coupling strength of the pendulum

        
    Attributes
    -----------
    times : 1d-array
        Stores time values for each time step.
    thetas : 1d-array
       Stores angular displacement values for each time step of the right pendulum.
    phis : 1d-array
       Stores angular displacement values for each time step of the left pendulum.
        
    Methods
    -----------
    kuramoto_model:
        Runs the Kuramoto model.
    kuramoto:
        The equation of the Kuramoto model.
    """

    def __init__(self, theta_0, phi_0, eta, n_iter, theta_length, phi_length, coupling_constant):
        self.theta_0 = theta_0
        self.phi_0 = phi_0
        self.eta = eta
        self.n_iter = n_iter
        self.theta_length = theta_length
        self.phi_length = phi_length
        self.coupling_constant = coupling_constant
        self.times = np.zeros(self.n_iter)
        self.thetas = np.zeros(self.n_iter)
        self.phis = np.zeros(self.n_iter)

    def kuramoto_model(self):
        self.thetas[0] = self.theta_0 * np.pi / 180.0
        self.phis[0] = self.phi_0 * np.pi / 180.0

        for i in range(self.n_iter - 1):
            self.times[i + 1] = self.times[i] + self.eta
            self.thetas[i + 1] = self.thetas[i] + self.eta * equation(self.thetas[i], self.phis[i], self.theta_length,
                                                                      self.coupling_constant)
            # print(self.equation(self.thetas[i], self.phis[i], self.theta_length, self.coupling_constant))
            self.phis[i + 1] = self.phis[i] + self.eta * equation(self.phis[i], self.thetas[i], self.phi_length,
                                                                  self.coupling_constant)
            # print(self.equation(self.phis[i], self.thetas[i], self.phi_length, self.coupling_constant))
        return self


system = Pendulum(theta_0=90, phi_0=30, eta=0.1, n_iter=100, theta_length=0.127, phi_length=0.127,
                  coupling_constant=0.5).kuramoto_model()
time = system.times * 180 / np.pi
right_pendulum = np.arcsin(np.sin(system.thetas)) * 180/np.pi
left_pendulum = np.arcsin(np.sin(system.phis)) * 180/np.pi

plt.plot(time, right_pendulum, lw=3, color='blue', label='theta')
plt.plot(time, left_pendulum, lw=3, color='red', label='phi')
plt.legend()
plt.xlabel('Time (s)', size=13)
plt.ylabel('Angle (deg)', size=13)
plt.title('Angle Against Time', size=13)
plt.show()
