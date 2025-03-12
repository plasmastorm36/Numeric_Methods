"""function comparing a bunch of numeric methods
"""
from typing import Callable
import numpy as np

def f(x: float)->float:
    """basic parabolic function

    Args:
        x (float): point on the x line

    Returns:
        float: value of y at x
    """
    return np.power(x, 2)

def midpoint (func: Callable[[float], float], a: float, b: float, n: int) -> float:
    """midpoint integration numeric method

    Args:
        func (callable[[float], float]): function
        a (float): begin
        b (float): end
        n (int): number of steps

    Returns:
        float: an approximation of the integral for any function from a to b
    """

    d = np.divide(np.subtract(b, a), n)
    s = 0.0
    x = np.divide(np.add(a, d), 2.0)
    while x < b:
        s = np.add(s, func(x))
        x = np.add(x, d)
    return np.multiply(s, d)

def trapezoid (func: Callable[[float], float], a: float, b: float, n: int) -> float:
    """trapezoid rule integration method

    Args:
        func (Callable[[float], float]): numeric function
        a (float): start
        b (float): end
        n (int): number of steps
        
    Returns:
        float: approximate integral of a function from a to b
    """

    d = np.divide(np.subtract(b, a), n)
    s = np.divide(np.add(func(a), func(b)), 2.0)
    x = np.add(a, d)
    while x < b:
        s = np.add(s, func(x))
        x = np.add(x, d)
    return np.multiply(s, d)

def simpson (func: Callable[[float], float], a: float, b: float, n: int) -> float:
    """simpson's rule integral approximation

    Args:
        func (Callable[[float], float]): any numerical function
        a (float): begin
        b (float): end
        n (int): steps

    Returns:
        float: approximation of definite integral of given function
    """
    if n < 2:
        print("n is too small")
        return
    if np.mod(n, 2) != 0:
        print("n is not even")
        return
    d = np.divide(np.subtract(b, a), n)
    s = np.add(func(a), func(b))
    x = np.add(a, d)
    i = 0
    while x < b:
        s = np.add(s, np.multiply(4.0, func(x))) if i % 2 == 0 else np.add(s,
        np.multiply(2.0, func(x)))
        x = np.add(x, d)
        i = np.add(1, i)
    return np.multiply(np.divide(d, 3), s)

if __name__ == '__main__':
    print(f"Midpoint: {midpoint(f, 0.0, 5.0, 10)}")
    print(f"Trapezoid: {trapezoid(f, 0.0, 5.0, 10)}")
    print(f"Simpson: {simpson(f, 0.0, 5.0, 10)}")
