import matplotlib.pyplot as plt
import numpy as np

# The two lines of data provided by the user
data_line1 = (
    "0.003346425 0.0  0.018800   0.000040   0.000248   0.000002   0.004325  0.000553"
)
exact_val1 = 0.0248812020
data_line2 = (
    "0.000022699 0.0  0.003759   0.000013   0.000008   0.000000   0.000349  0.000089"
)

# Create a dummy data file named 'data.txt'
try:
    with open("data.txt", "w") as f:
        f.write(data_line1 + "\n")
        f.write(data_line2 + "\n")
    print("Created data.txt with the provided data.")
except Exception as e:
    print(f"Error creating data.txt: {e}")


# Function to parse a data line
def parse_data_line(line_str):
    """Parses a string line into values and errors."""
    try:
        data_list = [float(x) for x in line_str.split()]
        # Values are at even indices (0, 2, 4, ...)
        val = data_list[0::2]
        # Errors are at odd indices (1, 3, 5, ...)
        err = data_list[1::2]
        if len(val) == len(err):
            return val, err
        else:
            return None, None
    except ValueError:
        return None, None


# Function to create and save a plot
def create_plot(orders, val, err, title, filename, exact_val):
    """Creates and saves a convergence plot."""
    # Create a new, separate figure for this plot
    plt.figure(figsize=(10, 6))

    # Plot the data with error bars
    plt.errorbar(
        orders, val, yerr=err, marker="o", capsize=5, linestyle="-", label=title
    )

    plt.axhline(y=exact_val, color="r", linestyle="--", label="Exact Value")

    # Add labels, title, and legend
    plt.xlabel("Order")
    plt.ylabel("Double Occupancy Value")
    plt.title(f"Convergence of Double Occupancy vs. Order\n({title})")
    plt.legend()

    # Set x-ticks to be integers for the orders
    if len(orders) > 0:
        plt.xticks(orders)

    # Add grid and ensure clean layout
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.tight_layout()

    # Save the plot to its own file
    plt.savefig(filename)
    print(f"Plot saved to {filename}")


# Main script to read file and generate plots
try:
    with open("double.txt", "r") as f:
        lines = f.readlines()

    if len(lines) < 2:
        print("Error: The file does not contain two lines of data.")
    else:
        # Process Dataset 1
        val1, err1 = parse_data_line(lines[0])
        if val1:
            orders1 = np.arange(len(val1))
            create_plot(
                orders1, val1, err1, "Dataset 1", "dataset1_convergence.png", exact_val1
            )
        else:
            print("Error parsing Dataset 1.")

        # Process Dataset 2
        val2, err2 = parse_data_line(lines[1])
        if val2:
            orders2 = np.arange(len(val2))
            create_plot(
                orders2, val2, err2, "Dataset 2", "dataset2_convergence.png", exact_val2
            )
        else:
            print("Error parsing Dataset 2.")

except FileNotFoundError:
    print("Error: data.txt not found.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
