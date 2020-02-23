import numpy as np
import matplotlib.pyplot as plt

def approx(N, r, b, ts):
    return np.exp(-r*N*np.exp(-r*ts) / b)
def actual(N, r, b, ts):
    return (b/r) * np.power(1 - r/(b*(np.exp(r*ts)-1)+r), N) *\
           r*np.exp(r*ts) / (b*(np.exp(r*ts)-1)+r)

# N, r, b = 1e8, 0.006554615384615381, 0.095585
# ts = np.linspace(0.6*np.log(N)/r, 1.3*np.log(N)/r, num = 1e4)
# R, = plt.plot(ts[:-1], np.diff(actual(N, r, b, ts)) / (ts[1]-ts[0]), label="line1")
# A, = plt.plot(ts[:-1], np.diff(approx(N, r, b, ts)) / (ts[1]-ts[0]), label="line2")
# plt.legend([R, A], ['Real','Approximate'])
# plt.xlabel("Time")
# plt.ylabel("Density")
# plt.title(r'Approximation level for stop time $\tau$')
# plt.show()

def testApprox(N, r, b):
    ts = np.linspace(0, 1.5*np.log(N)/r, num = 1e4)
    R, = plt.plot(ts[:-1], np.diff(actual(N, r, b, ts)) / (ts[1]-ts[0]), label="line1")
    A, = plt.plot(ts[:-1], np.diff(approx(N, r, b, ts)) / (ts[1]-ts[0]), label="line2")
    plt.legend([R, A], ['Real','Approximate'])
    plt.xlabel("Time")
    plt.ylabel("Density")
    plt.title(r'Approximation level for stop time $\tau$')
    plt.show()