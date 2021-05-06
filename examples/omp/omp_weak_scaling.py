import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import sys
import time

def main():
    threads = np.arange(1,9)
    for th in threads:
        print("Running experiment on {} thread(s) ...".format(th), end='')
        os.environ["OMP_NUM_THREADS"] = str(th)
        start = time.time()
        subprocess.run("./time_loop.exe")
        end = time.time()


if __name__ == '__main__':
    main()
