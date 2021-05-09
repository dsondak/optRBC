import sys
import matplotlib.pyplot as plt
import numpy as np

points = None
threads = []
times = []

f = open("multithread.txt", "r")

for line in f: 
    line = line.rstrip()
    if (line.startswith('N_exp')): 
        # Plot the dataset just finished
        if points is not None: 
            plt.plot(threads, times, label='2^' + str(points))
        threads = []
        times = []
        fields = line.split()
        points = fields[1]
    elif (line.startswith('Threads')): 
        thread_time = line.split(',')
        threads.append(int(thread_time[0].split()[-1]))
        times.append(float(thread_time[1].split()[-2]))

plt.legend()
plt.xlabel("Num Threads")
plt.ylabel("Execution time (s)")
plt.show()
f.close()