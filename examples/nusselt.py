import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

nu1 = pd.read_csv('nu1.txt', header=None)
nu2 = pd.read_csv('nu2.txt', header=None)
nu3 = pd.read_csv('nu3.txt', header=None)

time = np.arange(0,100,0.1)

f, ax = plt.subplots(figsize=(7,5))

ax.plot(time, nu1.values, label=r'Ra$=5\times 10^4$', linewidth=3, alpha=0.5)
ax.plot(time, nu2.values, label=r'Ra$=5\times 10^3$', linewidth=3, alpha=0.5)
ax.plot(time, nu3.values, label=r'Ra$=5\times 10^2$', linewidth=3, alpha=0.5)

ax.set_xlabel('Time (s)', size=15)
ax.set_ylabel('Nusselt number', size=15)
ax.legend()
plt.savefig('nusselt.png')