#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 22:54:04 2023

@author: jeremie
"""



import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt

# Define the LinearTransform module
class LinearTransform(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(LinearTransform, self).__init__()
        self.A = nn.Parameter(torch.randn(input_dim, output_dim)/6)

    def forward(self, x):
        projected_points = torch.matmul(x, self.A)
        return projected_points
    
    def compute_cost(self, x):
        projected_points = self.forward(x)
        norms = torch.norm(projected_points, dim=1)
        normalized_projected_points = projected_points / norms.unsqueeze(1)
        cost = torch.norm(projected_points - normalized_projected_points)
        return cost

def plot_projection(projected_points, time=None, title="Linear projection in 2D"):
    # Set the seed for reproducibility
    np.random.seed(42)

    # Check if time is provided
    if time is not None and len(time) != projected_points.shape[0]:
        raise ValueError("The length of time should match the number of projected points.")

    # Convert the projected points array to a tensor
    projected_points_tensor = torch.Tensor(projected_points)

    # Compute the norm of each projected point
    norms = torch.norm(projected_points_tensor, dim=1)

    # Normalize the projected points
    normalized_projected_points_tensor = projected_points_tensor / norms.unsqueeze(1)

    # Convert the normalized projected points tensor back to a numpy array
    normalized_projected_points = normalized_projected_points_tensor.detach().numpy()

    # Plot the projected points using scatter
    if time is not None:
        plt.scatter(projected_points[:, 0], projected_points[:, 1], label='Projected Points', s=5, c=time, cmap='viridis')
        plt.colorbar(label='Time')
    else:
        plt.scatter(projected_points[:, 0], projected_points[:, 1], label='Projected Points', s=5, c='blue')

    # Plot the normalized projected points using scatter
    plt.scatter(normalized_projected_points[:, 0], normalized_projected_points[:, 1],
                label='Normalized Projected Points', s=5, c='red')

    # Plot lines between projected points and normalized projected points
    for i in range(projected_points.shape[0]):
        plt.plot([projected_points[i, 0], normalized_projected_points[i, 0]],
                 [projected_points[i, 1], normalized_projected_points[i, 1]],
                 color='gray', alpha=0.3)

    # Set x and y axes to have the same scale
    plt.axis('equal')

    # Plot a circle with center (0, 0) and radius 1
    theta = np.linspace(0, 2 * np.pi, 100)
    x = np.cos(theta)
    y = np.sin(theta)
    plt.plot(x, y, 'r--')

    plt.title(title)
    plt.xlabel("Dimension 1")
    plt.ylabel("Dimension 2")
    plt.legend()



def plot_projection_before_after_optimisation(initial_projection, final_projection, time=None):
    # Set the seed for reproducibility
    np.random.seed(42)

    # Check if time is provided
    if time is not None and len(time) != initial_projection.shape[0]:
        raise ValueError("The length of time should match the number of projected points.")

    # Get the x and y limits for both initial and final projections
    x_limits = [min(np.min(initial_projection[:, 0]), np.min(final_projection[:, 0])),
                max(np.max(initial_projection[:, 0]), np.max(final_projection[:, 0]))]
    y_limits = [min(np.min(initial_projection[:, 1]), np.min(final_projection[:, 1])),
                max(np.max(initial_projection[:, 1]), np.max(final_projection[:, 1]))]

    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plot_projection(initial_projection, time, title="Linear transformation before optimisation")
    plt.xlim(x_limits)
    plt.ylim(y_limits)

    # Plot the final projection
    plt.subplot(1, 2, 2)
    plot_projection(final_projection, time, title="Linear transformation after optimisation")

    # Set the same x and y limits for both subplots
    plt.xlim(x_limits)
    plt.ylim(y_limits)

    plt.tight_layout()
    plt.show()

# Perform the linear transformation using the forward method
def optimize_linear_transform(points, time=None):
    # Convert points to a torch tensor
    points_tensor = torch.Tensor(points)

    # Create an instance of the LinearTransform module
    input_dim = points.shape[1]
    output_dim = 2
    linear_transform = LinearTransform(input_dim, output_dim)

    # Perform the linear transformation using the forward method
    initial_projection = linear_transform(points_tensor).detach().numpy()
    plot_projection(initial_projection, time)
    plt.show()

    # Compute the cost
    cost = linear_transform.compute_cost(points_tensor)
    print("Cost:", cost.item())

    # Perform automatic differentiation and compute gradients
    cost.backward()
    gradients = linear_transform.A.grad
    print("Gradients:")
    print(gradients)

    # Define the optimizer and learning rate
    optimizer = torch.optim.Adam(linear_transform.parameters(), lr=0.1)

    # Perform gradient descent
    num_iterations = 1000
    cost_history = []

    for i in range(num_iterations):
        # Zero the gradients
        optimizer.zero_grad()

        # Compute the cost
        cost = linear_transform.compute_cost(points_tensor)
        cost_history.append(cost.item())

        # Perform automatic differentiation
        cost.backward()

        # Update the parameters
        optimizer.step()

        # Print the cost every 100 iterations
        if (i + 1) % 100 == 0:
            print("Iteration", i + 1, "Cost:", cost.item())

    # Plot the cost evolution
    plt.plot(cost_history)
    plt.title("Cost Evolution")
    plt.xlabel("Iteration")
    plt.ylabel("Cost")
    plt.show()

    # Perform the linear transformation using the forward method again
    final_projection = linear_transform(points_tensor).detach().numpy()
    plot_projection(final_projection, time)
    plt.show()
    plot_projection_before_after_optimisation(initial_projection, final_projection, time)
    # Access the updated matrix A
    updated_A = linear_transform.A

    print("Updated Matrix A:")
    print(updated_A)


# Generate 100 points in 10 dimensions
np.random.seed(42)  # Set the random seed for reproducibility
m = 100  # Number of points
n = 40  # Number of dimensions
points = np.random.randn(m, n)/2

optimize_linear_transform(points, time=None)


#
np.random.seed(42)
num_points = 100
time_start = 0
time_end = 0.3 * np.pi
num_dimensions = 4
period = 2 * np.pi
noise_std = 0.4
# Generate time points
time_points = np.linspace(time_start, time_end, num_points)

# Generate phases for each dimension
phases = np.random.uniform(low=0, high=2 * np.pi, size=num_dimensions)

# Generate points from the multidimensional sinusoidal function
points = []
for t in time_points:
    point = 2 * np.sin(period * t + phases) + np.random.normal(0, noise_std, size=num_dimensions)
    points.append(point)
points = np.array(points)
optimize_linear_transform(points, time=time_points)
