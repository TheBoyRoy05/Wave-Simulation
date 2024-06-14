import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
length = 1
speed = 1
gamma = 1
resolution = 100
num_steps = 2000
wave_width = 0.1
damping = 0.1

dx = length / resolution
dt = dx / speed

x = np.arange(0, length + dx, dx)
f = np.zeros((len(x), num_steps))

# Starting conditions
f[:, 0] = 0.5*np.sin(2 * np.pi * x) + 0.5*np.sin(3 * np.pi * x)
# f[:, 0] = np.exp(-(x - length/2)**2 / wave_width**2)
f[1:-1,1] = f[1:-1,0] + 0.5 * speed**2 * (f[2:,0] + f[:-2,0] - 2*f[1:-1,0]) * (dt/dx)**2 

fig, ax = plt.subplots()
ax.set(xlim=[0, length], ylim=[-1.2, 1.2])
line, = ax.plot(x, f[:, 0])

for i in range(2, num_steps):
    lagrangian = f[2:, i-1] + f[:-2, i-1] - 2*f[1:-1, i-1]
    f[1:-1, i] = (4*f[1:-1, i-1] - f[1:-1, i-2] * (2 - gamma*dt) + speed**2 * \
        lagrangian * (dt/dx)**2) / (2 + gamma*dt)

def update(frame):
    line.set_data(x, f[:, frame])
    ax.set_title(f"time: {frame*dt:.2f}")
    return line,

anim = FuncAnimation(fig=fig, func=update, frames=num_steps, interval=dt)
plt.show()