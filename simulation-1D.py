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
f = np.zeros((num_steps, len(x)))

# Starting conditions
f[0, :] = -0.5*np.sin(2 * np.pi * x) - 0.3*np.sin(3 * np.pi * x)
f[0, :] = np.exp(-(x - length/3)**2 / wave_width**2)
f[1, 1:-1] = f[0, 1:-1] + 0.5 * speed**2 * (f[0, 2:] + f[0, :-2] - 2*f[0, 1:-1]) * (dt/dx)**2

# Boundary conditions (Dirichlet, Neumann, Absorbing, Periodic)
left_boundary = "Absorbing"
right_boundary = "Periodic"

def get_boundary(type, side, i):
    neighbor = int(side + (-1)**side)
    if type == "Absorbing": 
        return f[i, side] + speed * (f[i, neighbor] - f[i, side]) * dt / dx
    elif type == "Neumann":
        return f[i+1, neighbor]
    elif type == "Periodic":
        return (f[i+1, 1] + f[i+1, -2]) / 2
    else: # Dirichlet
        return 0

for i in tqdm(range(1, num_steps-1)):
    lagrangian = f[i, 2:] + f[i, :-2] - 2*f[i, 1:-1]
    f[i+1, 1:-1] = (4*f[i, 1:-1] - f[i-1, 1:-1] * (2 - gamma*dt) + speed**2 * lagrangian * (dt/dx)**2) / (2 + gamma*dt)
    f[i+1, 0] = get_boundary(left_boundary, 0, i)
    f[i+1, -1] = get_boundary(right_boundary, -1, i)

# Plotting + Animation
fig, ax = plt.subplots()
ax.set(xlim=[0, length], ylim=[-1.2, 1.2])
line, = ax.plot(x, f[0, :])

def update(frame):
    line.set_data(x, f[frame, :])
    ax.set_title(f"time: {frame*dt:.2f}")
    return line,

anim = FuncAnimation(fig=fig, func=update, frames=num_steps, interval=dt)
# anim.save("Images/Wave-1D.gif", writer=PillowWriter(fps=30))
plt.show()