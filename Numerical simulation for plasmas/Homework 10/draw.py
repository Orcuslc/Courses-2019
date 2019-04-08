import matplotlib.pyplot as plt

with open("pi_time.txt", "r") as f:
	data = f.read().split('\n')

time = [float(x) for x in data[1:-1:2]]
plt.semilogx([1, 2, 4, 8, 16], time, '*-')
plt.show()