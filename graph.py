import numpy as np
import matplotlib.pyplot as plt

def read_file(filePath):
  data = np.loadtxt(filePath)
  x = data[:, 0]  # 1列目
  y = data[:, 1]  # 2列目
  return x, y


filePath = ["data/si-si.dat","data/si-o.dat","data/o-o.dat"]

plt.figure(figsize=(10,6))
x_data = []
y_data = []

for path in filePath:
  x ,y = read_file(path)
  x_data.append(x)
  y_data.append(y)
  
plt.plot(x_data[0],y_data[0],label=f"Si-Si",color="red")
plt.plot(x_data[1],y_data[1],label=f"Si-O",color="green")
plt.plot(x_data[2],y_data[2],label=f"O-O",color="blue")



plt.title("Visualization of RDF")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.legend()
plt.grid(True)

plt.show()
