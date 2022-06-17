import matplotlib as mat
import matplotlib.pyplot as plt
import numpy as np
plt.switch_backend("TkAgg")
# plt.switch_backend("Qt5Agg")

mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12

# data = np.loadtxt("g2w.txt")
data = np.loadtxt("g4.dat")
grid = data[:, 0]
atom = data[:, 1]
dimer = data[:, 2]
tetramer = data[:, 3]

plt.plot(grid, atom, label="Atom")
plt.plot(grid, dimer, label="Dimer")
plt.plot(grid, tetramer, label="2x2 Tetramer")
plt.legend()
plt.xlim([0.0, max(grid)+1.0e-1])
# plt.ylim([0.0, 0.1])
# plt.yscale("log")
# plt.xlabel("$\omega_n/t$")
# plt.ylabel("$G_2^c(i\omega_n, r_1, r_1)$")
plt.xlabel("$\\tau \\cdot t$")
plt.ylabel("$G_4^c$")
plt.savefig("g4r1r1.pdf")
plt.show()
