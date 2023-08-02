#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 22:54:04 2023

@author: jeremie
"""

from sklearn.decomposition import PCA
import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt

# Define the LinearTransform module
class LinearTransform(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(LinearTransform, self).__init__()
        self.A = nn.Parameter(torch.randn(input_dim, output_dim)/3)

    def forward(self, x):
        projected_points = torch.matmul(x, self.A)
        return projected_points
    
    def compute_cost(self, x):
        projected_points = self.forward(x)
        norms = torch.norm(projected_points, dim=1)
        normalized_projected_points = projected_points / norms.unsqueeze(1)
        cost = torch.norm(projected_points - normalized_projected_points)
        return cost

def plot_projection(projected_points, time=None, ax=None, limits=None, title=None):
    if ax is None:
        ax = plt.gca()

    if title is not None:
        ax.set_title(title)



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
        sc = ax.scatter(projected_points[:, 0], projected_points[:, 1], s=5, c=time, cmap='viridis')
        fig = ax.get_figure()
        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label('Time')
    else:
        ax.scatter(projected_points[:, 0], projected_points[:, 1], label='Projected Points', s=5, c='blue')


    # Plot the normalized projected points using scatter
    ax.scatter(normalized_projected_points[:, 0], normalized_projected_points[:, 1],
                label='Normalized Projected Points', s=5, c='red')

    # Plot lines between projected points and normalized projected points
    for i in range(projected_points.shape[0]):
        ax.plot([projected_points[i, 0], normalized_projected_points[i, 0]],
                 [projected_points[i, 1], normalized_projected_points[i, 1]],
                 color='gray', alpha=0.3)



    # Plot a circle with center (0, 0) and radius 1
    theta = np.linspace(0, 2 * np.pi, 100)
    x = np.cos(theta)
    y = np.sin(theta)
    ax.plot(x, y, 'r--')

    ax.set_xlabel("Dimension 1")
    ax.set_ylabel("Dimension 2")
    ax.legend()

    ax.axis('equal')
    if limits is not None:
        print(limits)
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.set_xbound(lower=limits[0], upper=limits[1])
        ax.set_ybound(lower=limits[0], upper=limits[1])
        ax.set_aspect('equal', adjustable='box')

def plot_projection_before_after_optimisation(initial_projection, final_projection, pca_projection, time=None):
    # Set the seed for reproducibility
    np.random.seed(42)

    # Check if time is provided
    if time is not None and len(time) != initial_projection.shape[0]:
        raise ValueError("The length of time should match the number of projected points.")

    # Get the x and y limits for both initial and final projections
    x_limits = [min(np.min(initial_projection[:, 0]), np.min(final_projection[:, 0]), np.min(pca_projection[:, 0])),
                max(np.max(initial_projection[:, 0]), np.max(final_projection[:, 0]), np.max(pca_projection[:, 0]))]
    y_limits = [min(np.min(initial_projection[:, 1]), np.min(final_projection[:, 1]), np.min(pca_projection[:, 1])),
                max(np.max(initial_projection[:, 1]), np.max(final_projection[:, 1]), np.max(pca_projection[:, 1]))]
    # limits = [1.1 * min(x_limits[0], y_limits[0]), 1.1 * max(x_limits[1], y_limits[1])]

    limits = [-1.1 * max(abs(x_limits[0]), abs(x_limits[1]), abs(y_limits[0]), abs(y_limits[1])),
              1.1 * max(abs(x_limits[0]), abs(x_limits[1]), abs(y_limits[0]), abs(y_limits[1]))]

    fig, axs = plt.subplots(1, 3, figsize=(20, 6))

    coeff_corr = np.sum(np.linalg.norm(final_projection, axis=1))/np.sum(np.linalg.norm(initial_projection, axis=1))
    axs[0].set_title("Linear transformation before optimization")
    sc = plot_projection(coeff_corr * initial_projection, time, ax=axs[0], limits=limits)

    axs[1].set_title("Linear transformation after optimization")
    sc = plot_projection(final_projection, time, ax=axs[1], limits=limits)

    coeff_corr = np.sum(np.linalg.norm(final_projection, axis=1))/np.sum(np.linalg.norm(pca_projection, axis=1))
    axs[2].set_title("PCA linear projection")
    sc = plot_projection(coeff_corr * pca_projection, time, ax=axs[2], limits=limits)

    fig.tight_layout()
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

    pca_projection = PCA(n_components=2).fit_transform(points)

    plot_projection_before_after_optimisation(initial_projection, final_projection, pca_projection, time)
    # Access the updated matrix A
    updated_A = linear_transform.A

    print("Updated Matrix A:")
    print(updated_A)


# Generate 100 points in 40 dimensions
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
num_dimensions = 30
period = 2 * np.pi
noise_std = 0.3
amp = 0.6
# Generate time points
time_points = np.linspace(time_start, time_end, num_points)

# Generate phases for each dimension
phases = np.random.uniform(low=0, high=2 * np.pi, size=num_dimensions)

# Generate points from the multidimensional sinusoidal function
points = []
for t in time_points:
    point = amp * np.sin(period * t + phases) + np.random.normal(0, noise_std, size=num_dimensions)
    points.append(point)
points = np.array(points)
optimize_linear_transform(points, time=time_points)
