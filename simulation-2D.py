import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation, PillowWriter
from tqdm import tqdm

# Constants
length = 1
speed = 0.1
gamma = 0
resolution = 100
num_steps = 1000
playback_speed = 2
CFL = np.sqrt(2)/2

dx = length / resolution
dt = CFL * dx / speed

x = np.arange(0, length + dx, dx)
y = np.arange(0, length + dx, dx)
X, Y = np.meshgrid(x, y)
f = np.zeros((num_steps, len(x), len(x)))

# Starting conditions
f[0, :, :] = np.exp(-((X - length/2)**2 + (Y - length/3)**2) / 0.05**2)
f[1, 1:-1, 1:-1] = f[0, 1:-1, 1:-1] + 0.5 * speed**2 * (dt/dx)**2 * \
    (f[0, 2:, 1:-1] + f[0, :-2, 1:-1] - 2*f[0, 1:-1, 1:-1] + 
     f[0, 1:-1, 2:] + f[0, 1:-1, :-2] - 2*f[0, 1:-1, 1:-1])

# Boundary conditions (Dirichlet, Neumann, Absorbing, Periodic)
top_boundary = "Dirichlet"
left_boundary = "Neumann"
right_boundary = "Absorbing"
bottom_boundary = "Periodic"

sides = {"top": "f[i, -1, :]", "left": "f[i, :, 0]", "right": "f[i, :, -1]", "bottom": "f[i, 0, :]"}

def get_neighbor(side, index="i"):
    side_info = sides[side][2:-1].split(', ')
    neighbor_index = lambda x: int(int(x) + (-1)**int(x))
    if side_info[1] == ':':
        return f"f[{index}, :, {neighbor_index(side_info[2])}]"
    else: 
        return f"f[{index}, {neighbor_index(side_info[1])}, :]"

def update_boundary(type, side, i, g=None):
    if type == "Absorbing": 
        return eval(sides[side]) + speed * (eval(get_neighbor(side)) - eval(sides[side])) * dt / dx
    elif type == "Neumann":
        return get_neighbor(side, "i+1") + dx * g(dt * i) if g else eval(get_neighbor(side, "i+1"))
    elif type == "Periodic":
        keys = list(sides.keys())
        opposite = get_neighbor(keys[3 - keys.index(side)], "i+1")
        return (eval(get_neighbor(side, "i+1")) + eval(opposite)) / 2
    else: # Dirichlet
        return 0

for i in tqdm(range(1, num_steps-1)):
    lagrangian = f[i, 2:, 1:-1] + f[i, :-2, 1:-1] + f[i, 1:-1, 2:] + f[i, 1:-1, :-2] - 4*f[i, 1:-1, 1:-1]
    f[i+1, 1:-1, 1:-1] = (4*f[i, 1:-1, 1:-1] - f[i-1, 1:-1, 1:-1] * (2 - gamma*dt) + \
        2 * speed**2 * (dt/dx)**2 * lagrangian) / (2 + gamma*dt)
    f[i+1, -1, :] = update_boundary(top_boundary, "top", i)
    f[i+1, :, 0] = update_boundary(left_boundary, "left", i)
    f[i+1, :, -1] = update_boundary(right_boundary, "right", i)
    f[i+1, 0, :] = update_boundary(bottom_boundary, "bottom", i)

# Plotting + Animation
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surface = ax.plot_surface(X, Y, f[0, :, :])
ax.set(xlim=[0, length], ylim=[0, length], zlim=[-0.5, 0.5])

def update(frame):
    global surface
    surface.remove()
    surface = ax.plot_surface(X, Y, f[int(frame*playback_speed), :, :], cmap=cm.coolwarm)
    ax.set_title(f"time: {frame*dt:.2f}")

anim = FuncAnimation(fig=fig, func=update, frames=int(num_steps//playback_speed), interval=dt*playback_speed)
# anim.save("Images/Wave-2D.gif", writer=PillowWriter(fps=30))
plt.show()