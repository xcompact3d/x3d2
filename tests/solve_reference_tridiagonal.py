import numpy as np

low = [1. / 4., 1. / 3., 1. / 3., 1. / 4., 2.]
up = [2., 1. / 4., 1. / 3., 1. / 3., 1. / 4.]

A = np.diag(low, k=-1) + np.diag(up, k=+1) + np.diag([1.] * 6)
b = np.array([0, 1, 3, 3, 2, 1])
X = np.linalg.solve(A, b)

a = 1. / 3.
A = np.diag([a] * 5, k=-1) + np.diag([a] * 5, k=+1) + np.diag([1.] * 6)
A[5, 0] = a
A[0, 5] = a
X = np.linalg.solve(A, b)
print(A)
