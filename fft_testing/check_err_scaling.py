import numpy as np
import matplotlib.pyplot as plt

c2c_errs = np.loadtxt('errs_c2c.txt')
r2c_errs = np.loadtxt('errs_r2c.txt')

fig, ax = plt.subplots(1,1, figsize=(10,6), constrained_layout=True)
size=22

save = True

ax.plot(c2c_errs[:,0], c2c_errs[:,1], ls='', marker='o', ms=12, label='C2C')
ax.plot(r2c_errs[:,0], r2c_errs[:,1], ls='', marker='^', ms=12, label='R2C')

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel(r'$N_{x}$', fontsize=size)
ax.set_ylabel(r'$e_{\infty}$', fontsize=size)

ax.tick_params(axis='both', labelsize=size)

ax.legend(fontsize=size)

if save:
   fig.savefig('comp_errs.png', dpi=300)

plt.show()
