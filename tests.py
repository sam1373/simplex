from simplex import Simplex
import numpy as np

def test(A, b, c):
    print("LP:")
    print("A:")
    print(A)
    print("b:")
    print(b)
    print("c:")
    print(c)

    sim = Simplex(A, b, c)

    feas = sim.isFeasible()
    print("This LP is", "not" * (not feas), "feasible.")

    if feas:
        sim.solve()
        print("Solution:")
        soln, cost = sim.getRes()
        print(soln)
        print("Objective cost:")
        print(cost)

#tests

print("-------------")

A = np.array([[1, 0, 0, 1, -2, 1], [0, 1, 0, 2, 0, -2], [0, 0, 1, -1, 1, 0]])
b = np.array([-2, 0, -1])
c = np.array([0, 0, 0, -2, -2, -1])

test(A, b, c)

print("-------------")

A = np.array([[1, 0, 0, -2, 3, 1], [0, 1, 0, 1, 3, -1], [0, 0, 1, 2, -1, -3]])
b = np.array([-1, -1, -6])
c = np.array([0, 0, 0, 1, 1, 1])

test(A, b, c)

print("-------------")

A = np.array([[1, 0, 0, 3, -2, 1], [0, 1, 0, 2, 1, -1], [0, 0, 1, 1, -1, 1]])
b = np.array([3, 2, 1])
c = np.array([0, 0, 0, 3, 2, 1])

test(A, b, c)

print("-------------")

A = np.array([[1, 0, 0, 2, 1, 0, 1], [0, 1, 0, -2, -1, 3, 4], [0, 0, 1, 2, -2, 1, 1]])
b = np.array([3, -2, -4])
c = np.array([0, 0, 0, -1, 1, 1, 4])

test(A, b, c)

print("-------------")

A = np.array([[1, 0, 0, 5, -4, -1], [0, 1, 0, 1, -1, 0], [0, 0, 1, -3, 4, 1]])
b = np.array([10, 4, 1])
c = np.array([1, 1, -1, 4, 9, 3])

test(A, b, c)

print("-------------")

A = np.array([[1, 0, 0, 1, 0, 1], [0, 1, 0, -1, -1, 1], [0, 0, 1, 1, 1, 1]])
b = np.array([2, 1, 3])
c = np.array([1, 1, -1, 3, 2, 1])

test(A, b, c)