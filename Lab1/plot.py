import numpy as np
import matplotlib.pyplot as plt

def f(x, const):
    return lambda x: An(n) * np.cos(const * n) + Bn(n) * np.sin(const * n)

def fab():
    return lambda n: 1/(np.pi**2 * n**2), lambda n: -1/(np.pi * n)

def plotSeries(x, S):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, S)
    plt.show()
   
a = 0
b = 1
a0 = 2/3
rng = int(input('Wpisz range: '))
T = b-a
x = np.linspace(a, b, num=1000)
S = np.zeros((x.shape))

An, Bn = fab()
const = 2 * np.pi * x / T
f = f(x, const)
import time

fig, ax = plt.subplots()
line, = ax.plot(x, S)
#ax.set_xlim(a,b)
#ax.set_ylim(-0.2,1.2)
ax.grid()
#plt.autoscale()
for k in range(1, rng+1):
    S = 0
    nrng = k
    for i in range(1, nrng+1):
        n = i
        S += f(x)
    S += a0/2
    line.set_ydata(S)
    ax.relim()
    ax.autoscale()
    plt.draw()
    plt.pause(0.05)


plt.show()

