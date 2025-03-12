"""basic file to handle numeric methods for derivatives
"""
import numpy as np
import matplotlib.pyplot as plt
def f(x: float) -> float:
    """ mathematic function to be used

    Args:
        x (float): variable to be calculated

    Returns:
        float: the function for x
    """
    return np.power(x, 2)

def backward(x: float, h: float) -> float:
    """back difference numerical method

    Args:
        x (float): point
        h (float): step size

    Returns:
        float: derivative at x
    """
    return np.divide(np.subtract(f(x), f(x - h)), h)

def forward(x: float, h: float) -> float:
    """forward difference numerical method

    Args:
        x (float): position
        h (float): step size

    Returns:
        float: derivative of function at point x
    """
    return np.divide(np.subtract(f(x + h), f(x)), h)

def central(x: float, h: float) -> float:
    """central difference numerical method

    Args:
        x (float): point
        h (float): step size

    Returns:
        float: derivative at point x
    """
    return np.divide(np.subtract(f(x + h), f(x - h)), np.multiply(2, h))

def richard(x: float, h: float) -> float:
    """extrapolation using other methods to provide more accuracy

    Args:
        x (float): point
        h (float): step size

    Returns:
        float: derivative at point x
    """
    return np.divide(np.subtract(np.multiply(4, central(x, h / 2)), central(x, h)), 3)

if __name__ == "__main__":
    points = np.arange(-5, 5, 0.1)
    solution = np.multiply(2, points)
    backward_solution = backward(points, 0.5)
    forward_solution = forward(points, 0.5)
    central_solution = central(points, 0.5)
    richard_solution = richard(points, 0.5)

    plt.plot(points, solution, color='blue', label='true derivative')
    plt.plot(points, backward_solution, color='purple', label='backward difference')
    plt.plot(points, forward_solution, color='orange', label = 'forward difference')
    plt.plot(points, central_solution, color='green', label = 'central difference')
    plt.plot(points, richard_solution, color ='red', label='richardson extrapolation')
    plt.title('numeric differentiation comparison')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.show()
