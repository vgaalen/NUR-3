import numpy as np

from ex1 import multid_quick_sort

# test multid_quick_sort
def f(x):
    return x

x = np.array( [[1],[4],[2]])

print(multid_quick_sort(f,x))