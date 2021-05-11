"""
heat_equation_ntd.py

@author: andrewbrettin
"""

import numpy as np
from numpy import array, arange, linspace, logspace, zeros, ones, append
from numpy import pi, e, exp, cos, sin, log, sqrt
from numpy.linalg import solve

from scipy.integrate import quad
import xarray as xr

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class HeatEquation:
    """
    Model for numerically solving the heat equation with von Neumann boundary conditions.
    .. math::
        u_t = u_xx, x \in [a,b]
        u(x,0) = u_0(x)
        u_x(a, t) = f_1(t)
        u_x(b, t) = f_2(t)

    The solution is expected to be compactly supported on [a,b], and the von Neumann 
    boundary conditions are handled via the Neumann to Dirichlet transform.

    Parameters:
        a : float
            Left endpoint of domain.
        b : float
            Right endpoint of domain.
        m : int
            Number of interior gridpoints.
        dx : float
            Grid spacing.
        T : float
            Total timespan of simulation
        dt : float
            Timestep.
        r : float
            Courant number 0.5 dt / dx^2.
        N : int
            Number of temporal points.
        xs : numpy.array
            List of spatial gridpoints, including the boundaries.
        ts : numpy.array
            List of temporal points.
        u0 : function
            Sets the initial data for u.
        ux_a : function
            Sets the Von Neumann boundary condition on the left endpoint.
        ux_b
            Sets the Von Neumann boundary condition on the right endpoint.
        u : xarray.DataArray
            2D grid containing the numerically computed values of u at points 
            xs, ts. The first index is for the spatial index and the second 
            index is for the temporal index.
        utrue : xarray.DataArray
            The true value of the solution

    """

    def __init__(self, a=0, b=1, m=999, T=0.5, dt=None, 
            u0=None, neumann=True, ux_a=lambda t : 0, ux_b=lambda t : 0):
        self.a = a
        self.b = b
        self.m = m
        self.dx = (self.b-self.a)/(self.m+1)
        self.T = T
        if dt == None:
            self.dt = self.dx
        else:
            self.dt == dt
        self.r = 0.5 * self.dt / self.dx**2
        self.N = int(self.T/self.dx)
        if u0 == None:
            t0 = -1/(16*log(np.finfo(float).eps))
            self.u0 = lambda x : exp(-(x - 1/2)**2/(4*t0))
        else:
            self.u0 = lambda x : 0
        self.xs = linspace(self.a, self.b, self.m+2)
        self.ts = arange(0, self.T, self.dt)

        self.neumann = neumann
        self.ux_a = ux_a
        self.ux_b = ux_b

        self.u = xr.DataArray(np.zeros((self.m+2, self.N)), dims=('x', 't'),
                coords={'x': self.xs, 't': self.ts})
        self.u[{'t':0}] = self.u0(self.xs)
    
    def __I_H(self, tk, u_x):
        '''
        Returns the history part of the integral I_H using the trapezoidal rule.
        Parameters:
            tk : np.array
                integration abscissas, ignoring the last element (which is handled by I_L).
            u_x : function
                Neumann condition as a function of time.
        '''
        f = lambda t : u_x(t)/sqrt(tk[-1] - t)
        if len(tk) <= 2:
            return 0
        else:
            dt = tk[1] - tk[0]
            return dt * (0.5*f(tk[0]) + 0.5*f(tk[-2]) + sum([f(t) for t in tk[1:-2]]))

    def __I_L(self, tk, u_x):
        '''
        Returns the local part of the integral, I_L, by approximating u_x as a linear function.
        '''
        if len(tk) == 1:
            return 0
        else:
            dt = tk[-1] - tk[-2]
            return sqrt(dt) * (4/3*u_x(tk[-1]) + 2/3*u_x(tk[-2]))

    def run(self):
        '''
        Runs the simulation and populates the grid u.
        Test 20:39 6 May 2021: Homogenous
        '''
        A1 = np.diag((1 + 2*self.r)*ones(self.m), 0) + np.diag(-self.r*ones(self.m-1), 1) + np.diag(-self.r*ones(self.m-1), -1)
        A2 = np.diag((1 - 2*self.r)*ones(self.m), 0) + np.diag(self.r*ones(self.m-1), 1) + np.diag(self.r*ones(self.m-1), -1)
        
        # Set Neumann to Dirichlet boundary conditions:
        if not self.neumann:
            # Set homogeneous Dirichlet boundary conditions
            u_a = lambda t : 0
            u_b = lambda t : 0
        else: 
            def u_a(t):
                tk = self.ts[self.ts <=t]
                return -1/sqrt(pi)*(self.__I_H(tk, self.ux_a) + self.__I_L(tk, self.ux_a))
            def u_b(t):
                tk = self.ts[self.ts <=t]
                return 1/sqrt(pi)*(self.__I_H(tk, self.ux_b) + self.__I_L(tk, self.ux_b))
        for j in range(0, self.N-1):
            u_last = np.array(self.u.isel(t=j))[1:-1] # Interior last values
            # Set RHS
            b = A2 @ u_last
            b[0] = b[0] + self.r * (u_a(self.ts[j]) + u_a(self.ts[j+1]))
            b[self.m-1] = b[self.m-1] + self.r * (u_b(self.ts[j]) + u_b(self.ts[j+1]))
            # Solve
            u_next = solve(A1, b)
            self.u[:,j+1] = [u_a(self.ts[j+1]), *u_next, u_b(self.ts[j+1])]

    def plot(self, ax=None, title=None):
        if ax == None:
            fig, ax = plt.subplots()
        self.u.isel(t=0).plot(ax=ax, label='Initial value')
        for j in [1, 5, 10, 50, 100, -1]:
            self.u.isel(t=j).plot(ax=ax, label='Numerical sol, t={}'.format(self.ts[j]))
        ax.legend(bbox_to_anchor=[1,0.5], loc='center left')
        ax.set_title(title)

    def return_true_solution(self )

    def plot_true_sol(self, N_max=42, ax=None, title=None):
        u = quad(self.u0,0,1)[0]*np.ones((len(self.xs), len(self.ts)))
        for n in range(1, N_max):
            f = lambda y : self.u0(y)*cos(n*pi*y)
            u = u + 2*quad(f,0,1)[0]*np.array([[cos(n*pi*x)*exp(-(n*pi)**2 *t) for t in self.ts] for x in self.xs])
        u = xr.DataArray(u, dims=('x', 't'), coords={'x':self.xs, 't':self.ts})
        if ax == None:
            fig, ax = plt.subplots()
        u.isel(t=0).plot(ax=ax, label='Initial value')
        for j in [1, 5, 10, 50, 100, -1]:
            u.isel(t=j).plot(ax=ax, label='True sol, t={}'.format(self.ts[j]))
        ax.legend(bbox_to_anchor=[1,0.5], loc='center left')
        ax.set_title(title)

    def saveanimation(self, filename, filetype='.gif'):
        '''
        Saves an animation of the simulation in gif form.
        Untested 20:44 6 May 2021
        '''
        plt.rcParams['animation.ffmpeg_path'] = '/Users/andrewbrettin/opt/anaconda3/bin/ffmpeg'
        fig, ax = plt.subplots()
        ax.set(xlim=(0,1), ylim=(0,1))
        data = self.u.data
        line = ax.plot(self.xs, data[:,0])[0]

        def animate(j):
            ax.set_xlabel('x')
            ax.set_ylabel('u')
            ax.set_title('time={}'.format(np.round(self.ts[j], 3)))
            line.set_ydata(data[:,j])

        anim = FuncAnimation(fig, animate, interval=25, frames=len(self.ts)-1, repeat=False)
        plt.draw()
        plt.show()
        anim.save(filename + filetype, writer="imagemagick")