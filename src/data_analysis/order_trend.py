import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import Polynomial

# Load data from v_order_analysis.csv sepearted by semicolon
data = pd.read_csv('v_order_analysis.csv', sep=';')

# Extract x and y values
x = data['particles']
y = data['order_parameter']

# Create a box plot
sns.boxplot(x='particles', y='order_parameter', data=data)

# Fit the data with a 3rd-degree polynomial
poly_fit = Polynomial.fit(x, y, 3)

# Add the fitted curve to the plot
x_fit = np.linspace(x.min(), x.max())
plt.plot(x_fit, poly_fit(x_fit), color='red', linewidth=2)

# Set labels for the axes
plt.ylabel('Order Parameter')
plt.xlabel('Particles')

# Display the plot
plt.show()