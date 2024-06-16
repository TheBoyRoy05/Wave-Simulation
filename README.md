# Wave Simulation
The wave equation is a second-order linear partial differential equation for the description of waves or standing wave fields such as mechanical waves (e.g. water waves, sound waves and seismic waves) or electromagnetic waves (including light waves). It arises in fields like acoustics, electromagnetism, and fluid dynamics.  

The scalar wave equation is

$$\frac{\partial^2 \psi}{\partial t^2} = c^2 \nabla^2 \psi$$

where $c$ is the speed of the wave and $\nabla^2$ is the laplacian operation (the divergence of the gradient). More specifically, $\nabla^2$ can be written as the following in $n$ dimensions.

$$\nabla^2 \psi = \sum_{i=0}^n \frac{\partial^2 \psi}{\partial x_i^2}$$

In this project, I will use the finite-difference method to solve the wave equation and create a simulation for it in both 1 and 2 dimensions using python with a variety of different boundary conditions. But first, here's a quick introduction to discretization.

## Discretization

In applied mathematics, discretization is the process of transferring continuous functions, models, variables, and equations into discrete counterparts. To discretize our partial derivatives, let's first consider how we would descritze $\psi$ with respect to one variable. To do this, we first need to split $\psi$ into finite pieces of width $\Delta x$. When we do this, we'll get something like:

<p align="center"><img src="Images/Discretization.jpg" alt="Discretization" width="400"/></p>

If we want to find $\displaystyle\frac{d \psi}{dx} \Big|_{\displaystyle x_i} = \psi'(x_i)$, we can calculate the average slope of the two lines next to $x_i$. Doing this, we find that 

$$\psi'(x_i) = \frac{1}{2}\left( \frac{\psi(x_{i+1}) - \psi(x_i)}{\Delta x} + \frac{\psi(x_i) - \psi(x_{i-1})}{\Delta x} \right) = \frac{\psi(x_{i+1}) - \psi(x_{i-1})}{2 \Delta x}$$

Now, if we want to find $\displaystyle\frac{d^2 \psi}{dx^2} \Big|_{\displaystyle x_i} = \psi''(x_i)$, we can calculate the change in slope on the two lines next to $x_i$. Doing this, we find that

$$\psi''(x_i) = \frac{1}{\Delta x}\left( \frac{\psi(x_{i+1}) - \psi(x_i)}{\Delta x} - \frac{\psi(x_i) - \psi(x_{i-1})}{\Delta x} \right) = \frac{\psi(x_{i+1}) - 2\psi(x_i) + \psi(x_{i-1})}{\Delta x^2}$$

These two equations will give us what we need to apply the finite-difference method to the wave equation. Now, we're ready to solve the one dimensional wave equation.

## 1D Wave Equation

In one dimension, our wave equation becomes 

$$\frac{\partial^2 \psi}{\partial t^2} = c^2 \frac{\partial^2 \psi}{\partial x^2}$$

where $\psi$ is a function of both $t$ and $x$. In order to descritize $\psi$, I will write it as $\psi(t_i, x_j)$ where $x_j$ is a spatial coordinate and $t_i$ is a temporal coordinate. Applying the finite-difference method, we can turn the wave equation into

$$\frac{\psi(t_{i+1}, x_j) - 2\psi(t_i, x_j) + \psi(t_{i-1}, x_j)}{\Delta t^2} = c^2\frac{\psi(t_i, x_{j+1}) - 2\psi(t_i, x_j) + \psi(t_i, x_{j-1})}{\Delta x^2}$$

Our goal is to solve for $\psi(t_{i+1}, x_j)$ since we want to know what the next state of a cell is given it's previous states and the states of the cells nearby. To do this, we can simply rearrange our equation to get

$$\psi(t_{i+1}, x_j) = 2\psi(t_i, x_j) - \psi(t_{i-1}, x_j) + c^2\frac{\Delta t^2}{\Delta x^2} (\psi(t_i, x_{j+1}) - 2\psi(t_i, x_j) + \psi(t_i, x_{j-1}))$$

Hence, we have solved the wave equation in one dimension using the finite-difference method. Here's one such solution for the 1D wave function with a starting condition of a superposition of sine waves.

<p align="center"><img src="Images/Wave-1D-Sine.gif" alt="1D Sine Wave" width="400"/></p>

One thing to note is that $C = c\displaystyle\frac{\Delta t}{\Delta x}$ term is call the Courant number and the Courant-Friedrichs-Lewy (CFL) condition requires that $C \le 1$ for numerical stability for this simulation.

Now let's solve the two dimensional wave equation using the finite difference method and applying the same ideas that we used here.

## 2D Wave Equation

In two dimensions, our wave equation becomes 

$$\frac{\partial^2 \psi}{\partial t^2} = c^2\left( \frac{\partial^2 \psi}{\partial x^2} + \frac{\partial^2 \psi}{\partial y^2} \right)$$

where $\psi$ is a function of $x$, $y$, and $t$. In order to save space since the equation becomes quite large, I will write $\psi(t_i, x_j, y_k)$ as $\psi_{i,j,k}$ where each subscript corresponds to their respective coordinate. Here, $x_j$ and $y_k$ are our spatial coordinates and $t_i$ is our temporal coordinate. Applying the finite-difference method, we can turn the wave equation into

$$\frac{\psi_{i+1,j,k} - 2\psi_{i,j,k} + \psi_{i-1,j,k}}{\Delta t^2} = c^2\left(\frac{\psi_{i,j+1,k} - 2\psi_{i,j,k} + \psi_{i,j-1,k}}{\Delta x^2} + \frac{\psi_{i,j,k+1} - 2\psi_{i,j,k} + \psi_{i,j,k-1}}{\Delta y^2} \right)$$

For this simulation, $\Delta x = \Delta y$ to keep things simple. Now let's use this substitution to turn our equation into 

$$\frac{\psi_{i+1,j,k} - 2\psi_{i,j,k} + \psi_{i-1,j,k}}{\Delta t^2} = c^2\frac{\psi_{i,j+1,k} + \psi_{i,j-1,k} + \psi_{i,j,k+1} + \psi_{i,j,k-1} - 4\psi_{i,j,k}}{\Delta x^2}$$

Our goal is still to solve for $\psi_{i+1,j,k}$ since we want to know what the next state of a cell is given it's previous states and the states of the cells nearby. To do this, we can simply rearrange our equation to get

$$\psi_{i+1,j,k} = 2\psi_{i,j,k} - \psi_{i-1,j,k} + c^2\frac{\Delta t^2}{\Delta x^2}(\psi_{i,j+1,k} + \psi_{i,j-1,k} + \psi_{i,j,k+1} + \psi_{i,j,k-1} - 4\psi_{i,j,k})$$

Hence, we have solved the wave equation in two dimensions using the finite-difference method. Here's one such solution for the 2D wave function with a starting condition of a 2D Gaussian.

<p align="center"><img src="Images/Wave-2D-Gaussian.gif" alt="2D Gaussian" width="400"/></p>

## Damping

Now, let's solve the wave equation with damping by modifying the wave equation to include a damping term. This equation is given by

$$\frac{\partial^2 \psi}{\partial t^2} + \gamma \frac{\partial \psi}{\partial t} = c^2 \nabla^2 \psi$$

where $\gamma$ is the damping coefficient. Let's start in one dimension. Using the notation from the last part and the tools given to us by the finite-difference method, we can write this equation as 

$$\frac{\psi_{i+1,j} - 2\psi_{i,j} + \psi_{i-1,j}}{\Delta t^2} + \gamma \frac{\psi_{i+1,j} - \psi_{i-1,j}}{2\Delta t} = c^2\frac{\psi_{i,j+1} - 2\psi_{i,j} + \psi_{i,j-1}}{\Delta x^2}$$

We can start simplifying this by combining the left side into one fraction. This gives us 

$$\frac{2\psi_{i+1,j} - 4\psi_{i,j} + 2\psi_{i-1,j} + \gamma\psi_{i+1,j} - \gamma\psi_{i-1,j}}{2\Delta t^2} = c^2\frac{\psi_{i,j+1} - 2\psi_{i,j} + \psi_{i,j-1}}{\Delta x^2}$$

Now, if we combine like terms on the left side, we get 

$$\frac{(2 + \gamma)\psi_{i+1,j} - 4\psi_{i,j} + (2 - \gamma)\psi_{i-1,j}}{2\Delta t^2} = c^2\frac{\psi_{i,j+1} - 2\psi_{i,j} + \psi_{i,j-1}}{\Delta x^2}$$

Our goal is still to solve for $\psi_{i+1,j}$ since we want to know what the next state of a cell is given it's previous state and the states of the cells nearby. To do this, we can simply rearrange our equation to get

$$\psi_{i+1,j} = \frac{1}{2 + \gamma} \left( 4\psi_{i,j} - (2 - \gamma)\psi_{i-1,j} + 2c^2\frac{\Delta t^2}{\Delta x^2}(\psi_{i,j+1} - 2\psi_{i,j} + \psi_{i,j-1}) \right)$$

Hence we have solved the 1D wave equation with damping. Here's the same example as before, but with $\gamma = 1$. Notice how the wave gets dampened over time. This is exactly the result we hoped for.

<p align="center"><img src="Images/Wave-1D-Sine-Damped.gif" alt="Damped Sine Wave" width="400"/></p>

One thing to note is that if $\gamma = 0$, the $\frac{1}{2}$ at the start of the equation will cancel with the terms in the inside, giving us the original solution to the 1D wave equation. This confirms that our we did our math correctly.

## Boundary Conditions

If you were paying attention really closely, you might have noticed that the finite difference method doesn't tell us what to do at the boundaries. Let's take a look at the equations to see why.

$$\psi'(x_i) = \frac{\psi(x_{i+1}) - \psi(x_{i-1})}{2 \Delta x}$$

$$\psi''(x_i) = \frac{\psi(x_{i+1}) - 2\psi(x_i) + \psi(x_{i-1})}{\Delta x^2}$$

As it turns out, both the first and second derivatives rely on the cells before and after the cell we want to update. But what if there's no cell before, like the case of the first cell, or no cell after, like the case of the last cell. Turns out that there are many solutions to this problem, each with their own use case. I will demonstrate four such boundary conditions that are often used to solve this issue.

### Dirichlet Boundary Condition

This boundary condition, named after the German mathmatician Johann Peter Gustav Lejeune Dirichlet, is the one that I have been displaying in my examples so far. This boundary condition can be written explicitly as 

$$\psi \Big|_{\partial M} = 0$$

where $M$ is a manifold (the general term for surface) and $\partial M$ is the boundary of $M$. In one dimension $\partial M$ is at $x = 0$ and $x = L$ where $L$ is the length of the domain. The Dirichlet Boundary Condition takes the far simplest approach in my opinion, by just setting the boundaries to be zero so we can forget about them entirely.

### Neumann Boundary Condition

This boundary condition, named after John von Neumann, a Hungarian-American mathematician and physicist, was the most difficult to wrap my head around since it can be considerabley generalized. This condition is used when the rate of change of the wave function is known at the boundary. The most general form of this condition is 

$$\frac{\partial \psi}{\partial n} \Big|_{\partial M} = f(t)$$

where $\displaystyle\frac{\partial}{\partial n}$ is the outward pointing normal from a general boundary. In one dimension we have

$$\frac{\partial}{\partial n} \Big|_{x=L} = \frac{\partial}{\partial x},\ \ \ \ \frac{\partial}{\partial n} \Big|_{x=0} = -\frac{\partial}{\partial x}$$

Furthermore, $f(t)$ is the rate of the change of the wave function at the boundary. The simplest form of this boundary condition is when $f(t) = 0$, which is what I have implemented. In this case, there's no inward or outward flow and so the boundary cells will tend towards their neighboring cells.

In one dimension, this can be achieved by just setting the boundary cells to be equal to it's neighboring cell (Figure 1), or you could have update the boundary cell to move towards its neighbor by cutting the distance between them in half with each time step (Figure 2). The both achieve similar results as seen below, but the second method produces small waves on the left boundary which may not be desirable.

<figure align="center"><img src="Images/Wave-1D-Sine-Neumann1.gif" alt="Nuemann Method 1" width="400"/><figcaption align="center">Figure 1</figcaption></figure>

<figure align="center"><img src="Images/Wave-1D-Sine-Neumann2.gif" alt="Nuemann Method 1" width="400"/><figcaption align="center">Figure 2</figcaption></figure>

In two dimensions, we can't set the boundary cell to any one of it's neighbors since it it has multiple, so we have to use the second method of having it tend towards some average of the neighbors. [This article](https://folk.ntnu.no/leifh/teaching/tkt4140/._main056.html#kap:722) discusses one such way to do this using ghost cells around the boundary. However, it should be noted that they are applying it to the heat equation rather than the wave equation which is similar, but first order in time. Regardless, the Neumann boundary condition can be applied the same way. 

### Absorbing Boundary Condition (ABC)

ABCs are designed to simulate an open boundary where waves can exit the simulation domain without reflecting back. 