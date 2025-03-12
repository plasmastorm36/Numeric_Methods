""" basic script to compare accuracy of the RK methods
"""
from typing import Callable
import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import ArrayLike

def rk1 (f: Callable[[float, ArrayLike], ArrayLike], y: ArrayLike, t0: float, h: float,
            n: int)-> dict:
    """first order Runge Kutta method

    Args:
        f (Callable[[float, ArrayLike], ArrayLike]): ODE system representing the ODE
        y (ArrayLike): IVPs
        t0 (float): initial point for IVPs
        h (float): step size
        n (int): number of steps

    Returns:
        dict: dictionary to hold the t and y values
    """
    t = t0
    y_calc = np.copy(y)
    points = {'t': np.empty(n), 'y': np.empty(n)}
    for i in range(n):
        y_calc = np.add(y_calc, np.multiply(h, f(t, y_calc)))
        t = np.add(t, h)
        points['t'][i] = t
        points['y'][i] = y_calc[0]
    return points

def rk2 (f: Callable[[float, ArrayLike], ArrayLike], y: ArrayLike, t0: float, h: float,
            n: int) -> dict:
    """Second order Runge Kutta method

    Args:
        f (Callable[[float, ArrayLike], ArrayLike]): ODE System representing the ODE
        y (ArrayLike): IVPs
        t0 (float): initial point for IVP
        h (float): step size
        n (int): number of steps

    Returns:
        dict: dictionary containing t and solved y values
    """
    t = t0
    y_calc = np.copy(y)
    points = {'t': np.empty(n), 'y': np.empty(n)}
    for i in range(n):
        k1 = np.multiply(h, f(t, y_calc))
        k2 = np.multiply(h, f(np.add(t, h), np.add(y_calc, k1)))
        y_calc = np.add(y_calc, np.multiply(0.5, np.add(k1, k2)))
        t = np.add(t, h)
        points['t'][i] = t
        points['y'][i] = y_calc[0]
    return points

def rk3 (f: Callable[[float, ArrayLike], ArrayLike], y: ArrayLike, t0: float, h: float,
            n: int) -> dict:
    """Third order Runge Kutta method

    Args:
        f (Callable[[float, ArrayLike], ArrayLike]): ODE system representing ODE
        y (ArrayLike): IVPs
        t0 (float): point where IVPs occur
        h (float): step size
        n (int): number of steps

    Returns:
        dict: _description_
    """
    t = t0
    y_calc = np.copy(y)
    points = {'t': np.empty(n), 'y': np.empty(n)}
    for i in range(n):
        k1 = np.multiply(h, f(t, y_calc))
        k2 = np.multiply(h, f(np.add(t, np.multiply(0.5, h)), np.add(y_calc, np.multiply(0.5, k1))))
        k3 = np.multiply(h, f(np.add(t, h), np.add(np.subtract(y_calc, k1), np.multiply(2.0, k2))))
        y_calc = np.add(y_calc, np.divide(np.add(k1, np.add(np.multiply(4, k2), k3)), 6.0))
        t = np.add(t, h)
        points['t'][i] = t
        points['y'][i] = y_calc[0]
    return points

def rk4 (f: Callable[[float, ArrayLike], ArrayLike], y: ArrayLike, t0: float, h: float,
         n: int) -> dict:
    """Fourth Order Runge Kutte method

    Args:
        f (Callable[[float, ArrayLike], ArrayLike]): ODE system representing ODE
        y (ArrayLike): IVPs
        t0 (float): point where IVPs occur
        h (float): step size
        n (int): number of steps

    Returns:
        dict: dictionary containing point t and calculated value y
    """
    t = t0
    y_calc = np.copy(y)
    points = {'t': np.empty(n), 'y': np.empty(n)}
    for i in range(n):
        k1 = np.multiply(h, f(t, y_calc))
        k2 = np.multiply(h, f(np.add(t, np.multiply(0.5, h)), np.add(y_calc, np.multiply(0.5, k1))))
        k3 = np.multiply(h, f(np.add(t, np.multiply(0.5, h)), np.add(y_calc, np.multiply(0.5, k2))))
        k4 = np.multiply(h, f(np.add(t, h), np.add(y_calc, k3)))
        y_calc = np.add(y_calc, np.divide(np.add(k1, np.add(np.multiply(2, k2),
                np.add(np.multiply(2, k3), k4))),6))
        t = np.add(t, h)
        points['t'][i] = t
        points['y'][i] = y_calc[0]
    return points

def ode (t: float, y: ArrayLike) -> ArrayLike:
    """ODE system for a simple harmonic oscillator where omega is 1

    Args:
        t (float): point in time
        y (ArrayLike): current values for all derivatives and solution

    Returns:
        ArrayLike: all values for the ODE system
    """
    dy = np.empty(len(y))
    dy[0] = y[1]
    dy[1] = -y[0]
    return dy

def exact_solution (t0: float, h: float, n: int) -> dict:
    """exact solution to ODE above

    Args:
        t0 (float): initial solution point
        h (float): step size
        n (int): number of steps

    Returns:
        dict: dictionary of array of points t, and corresponding values y
    """
    points = {'t': np.empty(n), 'y': np.empty(n)}
    t = t0
    for i in range(n):
        np.add(t, h)
        points['t'][i] = t
        points['y'][i] = np.cos(t)
    return points

if __name__ == '__main__':
    y = np.empty(2)
    y[0] = 1
    y[1] = 0
    exact_points = exact_solution(0, 0.1, 100)
    rk1_points = rk1(ode, y, 0.0, 0.1, 100)
    rk2_points = rk2(ode, y, 0.0, 0.1, 100)
    rk3_points = rk3(ode, y, 0.0, 0.1, 100)
    rk4_points = rk4(ode, y, 0.0, 0.1, 100)

    plt.plot(exact_points['t'], exact_points['y'], color='green', label='exact')
    plt.plot(rk1_points['t'], rk1_points['y'], color='blue', label='rk1')
    plt.plot(rk2_points['t'], rk2_points['y'], color='purple', label='rk2')
    plt.plot(rk3_points['t'], rk3_points['y'], color='red', label='rk3')
    plt.plot(rk4_points['t'], rk4_points['y'], color='orange', label='rk4')
    plt.xlabel('t')
    plt.ylabel('y')
    plt.xticks()
    plt.yticks()
    plt.legend()
    plt.show()
    