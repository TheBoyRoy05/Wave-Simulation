import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

# Constants
length = 1
speed = 0.1
resolution = 100
num_steps = 1000
wave_width = 0.05

dx = length / resolution
dt = np.sqrt(2)/2 * dx / speed

x = np.arange(0, length + dx, dx)
y = np.arange(0, length + dx, dx)
X, Y = np.meshgrid(x, y)
f = np.zeros((len(x), len(x), num_steps))

# Starting conditions
f[:, :, 0] = np.exp(-((X - length/2)**2 + (Y - length/2)**2) / wave_width**2)
f[1:-1, 1:-1, 1] = f[1:-1, 1:-1, 0] + 0.5 * speed**2 * (dt/dx)**2 * \
    (f[2:, 1:-1, 0] + f[:-2, 1:-1, 0] - 2*f[1:-1, 1:-1, 0] + 
     f[1:-1, 2:, 0] + f[1:-1, :-2, 0] - 2*f[1:-1, 1:-1, 0])

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surface = ax.plot_surface(X, Y, f[:, :, 0])
ax.set(xlim=[0, length], ylim=[0, length], zlim=[-1, 1])

for i in tqdm(range(2, num_steps)):
    spacial_info = f[2:, 1:-1, i-1] + f[:-2, 1:-1, i-1] + f[1:-1, 2:, i-1] + f[1:-1, :-2, i-1] - 4*f[1:-1, 1:-1, i-1]
    f[1:-1, 1:-1, i] = 2*f[1:-1, 1:-1, i-1] - f[1:-1, 1:-1, i-2] + speed**2 * (dt/dx)**2 * spacial_info

def update(frame):
    global surface
    surface.remove()
    surface = ax.plot_surface(X, Y, f[:, :, frame], cmap=cm.coolwarm)
    ax.set_title(f"time: {frame*dt:.2f}")

anim = FuncAnimation(fig=fig, func=update, frames=num_steps, interval=dt)
plt.show()