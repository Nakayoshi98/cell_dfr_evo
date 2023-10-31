import itertools
import numpy as np
import sys
args = sys.argv

genes = int(args[1])
A = [0, 1]
a = np.array(list(itertools.product(A, repeat=genes)))[1:]
G = 1e4

np.random.seed(int(args[3]))
delta = (float(args[4]))


def f3(dim):
    while True:
        x = np.random.randn(dim)
        r = np.linalg.norm(x)
        if r != 0.0:
            return x / r


num_noise_LP = int(args[2])
random_noise = np.array([delta * abs(f3(genes))
                         for i in range(num_noise_LP - 2**genes + 1)])
a_e = np.array([[delta * a[k][i] / np.sqrt(sum(a[k]))
                 for i in range(genes)] for k in range(len(a))])
noise_set = np.concatenate([random_noise, a_e])
np.savetxt(
    "../perturbation_lib/N_noise_genes" + str(
        genes) +
    "_delta_" +
    args[4] +
    "_seed_" +
    args[3] +
    "_noise_set_N=" +
    args[2] +
    "_including_e.txt",
    noise_set)
