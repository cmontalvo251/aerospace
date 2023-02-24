##CENTER OF MASS CODE CREATED BY CHAT GPT - 2/24/2023

"""
In this program, we first define the positions of each component as a 2D NumPy array positions. We also define the weight of each component as a 1D NumPy array weights.

We then calculate the total weight of the aircraft by using the np.sum() function on the weights array. We also calculate the weight moments by multiplying each component's position by its weight and then summing across all components. This is done using the element-wise multiplication operator (*) and the np.sum() function with the axis parameter set to 0.

Finally, we calculate the center of mass by dividing the weight moments by the total weight. We print the result using an f-string.

Note that this program assumes that the position and weight arrays are of the same length, and that the position of each component is given as a 3D vector. If your data is in a different format, you may need to modify the code accordingly.
"""

import numpy as np

# Positions of components
positions = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

# Weights of components
weights = np.array([1000, 2000, 3000])

# Calculate total weight
total_weight = np.sum(weights)

# Calculate weight moments
weight_moments = np.sum(positions * weights.reshape(-1, 1), axis=0)

# Calculate center of mass
center_of_mass = weight_moments / total_weight

# Print result
print(f"Center of mass: {center_of_mass}")
