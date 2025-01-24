import numpy as np
import matplotlib.pyplot as plt

file_path = "rdf.dat"

data = np.loadtxt(file_path)

x = data[:, 0]  # 1列目
y = data[:, 1]  # 2列目

plt.figure(figsize=(8, 6))
plt.plot(x, y, label="Data", marker="o")

plt.title("Data Visualization")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.legend()
plt.grid(True)

plt.show()
