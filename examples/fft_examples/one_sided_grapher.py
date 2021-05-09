import sys
import matplotlib.pyplot as plt
import numpy as np

points = []
one_time = np.array([])
two_time = np.array([])
average = None
speedups = []

f = open("one_sided.txt", "r")

for line in f: 
	line = line.rstrip()
	if (line.startswith('N_points_exp')): 
		fields = line.split()
		curr_points = int(fields[-1])
		if not points: 
			points.append(curr_points)
		elif points[-1] is not curr_points: 
			points.append(curr_points)
			curr_speedups = two_time / one_time
			average = np.mean(curr_speedups)
			speedups.append(average)
		one_time = np.array([])
		two_time = np.array([])
	elif (line.startswith('1')): 
		curr_one_time = np.array([float(line.split()[-2])])
		one_time = np.append(one_time, curr_one_time)
	elif (line.startswith('2')): 
		curr_two_time = np.array([float(line.split()[-2])])
		two_time = np.append(two_time, curr_two_time)

curr_speedups = two_time / one_time
average = np.mean(curr_speedups)
speedups.append(average)

print(speedups)
print(points)

plt.plot(points, speedups)
plt.xlabel("2^N")
plt.ylabel("Speedup relative to two-sided FFT")
plt.legend()
plt.show()
f.close()