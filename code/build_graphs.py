import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
sns.set()
path = '../data/results/'

graph = []
for file in os.listdir(path):
    result = (np.load(os.path.join(path,file))['result']).tolist()
    graph.append((int(file.replace('ResultsmatrixVensembleSplit', '').replace('ResultsmatrixVprec', '').replace('.npz', '')), result))

graph = list(sorted(graph, key=lambda x:x[0]))

sns.lineplot(*zip(*graph[:-1]), marker="o")
plt.scatter(graph[-1][0],graph[-1][1],c='r',s=30, marker='x')
plt.xlabel('Reduced Matrix Size', fontsize=14)
plt.ylabel('MSE', fontsize=12)
plt.gca().legend(('Ensemble','tsvd'),loc='upper left')
plt.ylim(0.115, 0.2)
plt.savefig('../data/results/MSE_graph.png')
plt.show()