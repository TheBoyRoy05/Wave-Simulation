import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.animation import FuncAnimation, PillowWriter

# Constants
length = 1
speed = 1
gamma = 0
resolution = 100
num_steps = 1000
wave_width = 0.1
CFL = 1

dx = length / resolution
dt =  CFL * dx / speed

x = np.arange(0, length + dx, dx)
f = np.zeros((len(x), num_steps))

# Starting conditions
f[:, 0] = 0.5*np.sin(2 * np.pi * x) + 0.5*np.sin(3 * np.pi * x)
# f[:, 0] = np.exp(-(x - length/3)**2 / wave_width**2)
f[1:-1,1] = f[1:-1,0] + 0.5 * speed**2 * (f[2:,0] + f[:-2,0] - 2*f[1:-1,0]) * (dt/dx)**2

# Boundary conditions (Neumann, Dirichlet, Absorbing, Periodic)
boundary = "Periodic"

fig, ax = plt.subplots()
ax.set(xlim=[0, length], ylim=[-1.2, 1.2])
line, = ax.plot(x, f[:, 0])

for i in tqdm(range(2, num_steps)):
    lagrangian = f[2:, i-1] + f[:-2, i-1] - 2*f[1:-1, i-1]
    f[1:-1, i] = (4*f[1:-1, i-1] - f[1:-1, i-2] * (2 - gamma*dt) + speed**2 * lagrangian * (dt/dx)**2) / (2 + gamma*dt)

    if boundary == "Absorbing":
        f[0, i] = f[0, i-1] + speed * (f[1, i-1] - f[0, i-1]) * dt / dx
        f[-1, i] = f[-1, i-1] - speed * (f[-1, i-1] - f[-2, i-1]) * dt / dx
        continue
    elif boundary == "Neumann":
        # start = 2 * (f[1, i-1] - f[0, i-1]) # Method 2, makes small vibrations
        # end = 2 * (f[-2, i-1] - f[-1, i-1])
        f[0, i] = f[1, i]
        f[-1, i] = f[-2, i]
        continue
    elif boundary == "Periodic":
        start = f[1, i-1] + f[-1, i-1] - 2*f[0, i-1]
        end = f[0, i-1] + f[-2, i-1] - 2*f[-1, i-1]
    else:
        start = 0
        end = 0

    f[0, i] = (4*f[0, i-1] - f[0, i-2] * (2 - gamma*dt) + 2 * speed**2 * start * (dt/dx)**2) / (2 + gamma*dt)
    f[-1, i] = (4*f[-1, i-1] - f[-1, i-2] * (2 - gamma*dt) + 2 * speed**2 * end * (dt/dx)**2) / (2 + gamma*dt)

def update(frame):
    line.set_data(x, f[:, frame])
    ax.set_title(f"time: {frame*dt:.2f}")
    return line,

anim = FuncAnimation(fig=fig, func=update, frames=num_steps, interval=dt)
# anim.save("Images/Wave-1D.gif", writer=PillowWriter(fps=30))
plt.show()