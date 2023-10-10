import matplotlib.pyplot as plt

# Load the data from a file or create a list with the data if it's not in a file
data = [
    [1.91039, 1],
    [0.04109, 0.5],
    [0.0294474, 0.2],
    [0.00586565, 0.125],
    [0.00297486, 0.1]
]

# Separate the x and y data into separate lists
x = [row[1] for row in data]
y = [row[0] for row in data]

# Create a line plot
plt.plot(x, y)

# Set the plot title and axis labels
plt.title("Approximation of Error")
plt.ylabel("Error")
plt.xlabel("Dt")

# Set the y axis to a logarithmic scale
plt.yscale('log')


# Show the plot
plt.show()
