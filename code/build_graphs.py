import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
sns.set()

# path = '../data/results/results/ensemble/'
#
#
# for type in os.listdir(path):
#     graph = []
#     time_graph= []
#     type_path = path + type
#     for file in os.listdir(type_path):
#         if file[-3:] == 'npz':
#             result = (np.load(os.path.join(type_path, file))['result']).tolist()
#             time = (np.load(os.path.join(type_path, file))['time']).tolist()
#             graph.append((int(file.replace('ResultsmatrixVensembleSplit', '').replace('ResultsmatrixVprec', '').replace('.npz', '').replace(type,'')), result))
#             time_graph.append((int(file.replace('ResultsmatrixVensembleSplit', '').replace('ResultsmatrixVprec', '').replace('.npz', '').replace(type,'')), time))
#     print(graph)
#     graph = list(sorted(graph, key=lambda x:x[0]))
#     time_graph = list(sorted(time_graph, key=lambda x: x[0]))
#
#     sns.lineplot(*zip(*graph[:-1]), marker="o")
#     plt.scatter(graph[-1][0],graph[-1][1],c='r',s=30, marker='x')
#     plt.xlabel('Reduced Matrix Size', fontsize=14)
#     plt.ylabel('MSE', fontsize=12)
#     plt.title(type)
#     plt.gca().legend(('Ensemble','tsvd'),loc='upper left')
#     # plt.ylim(0.115, 0.2)
#     if type=='state':
#         plt.title('pollutant')
#     plt.savefig('../data/results/MSE_graph'+type+'.png')
#     plt.show()
#     plt.clf()
#     sns.lineplot(*zip(*time_graph[:-1]), marker="o")
#     plt.scatter(time_graph[-1][0], time_graph[-1][1], c='r', s=30, marker='x')
#     plt.xlabel('Reduced Matrix Size', fontsize=14)
#     plt.ylabel('seconds', fontsize=12)
#     plt.title(type)
#     plt.gca().legend(('Ensemble', 'tsvd'), loc='upper left')
#     # plt.ylim(0.115, 0.2)
#     if type=='state':
#         plt.title('pollutant')
#     plt.savefig('../data/results/time_graph'+type+'.png')
#     plt.show()
#     plt.clf()
#
# explained_variance_fp = '../data/results/results/explained_variance_Ch.npy'
# horizontal = np.cumsum(list(sorted(np.load(explained_variance_fp), reverse=True)))
# plt.xlabel('Number of components', fontsize=14)
# plt.ylabel('Explained Variance', fontsize=12)
# plt.title('Horizontal Localisation')
# sns.lineplot(data=horizontal)
# plt.savefig('../data/results/explained_h_variance.png')
# plt.show()
# plt.clf()
# explained_variance_fp = '../data/results/results/explained_variance_Cv.npy'
# vertical = np.cumsum(list(sorted(np.load(explained_variance_fp), reverse=True)))
# plt.xlabel('Number of components', fontsize=14)
# plt.ylabel('Explained Variance', fontsize=12)
# plt.title('Vertical Localisation')
# sns.lineplot(data=vertical)
# plt.savefig('../data/results/explained_v_variance.png')
# plt.show()
# plt.clf()

# msexB = np.load('../data/results/results/distributed_in_time_rh3_200matrixVensembleSplit40state.npz')['msexB']
# print(msexB)
# msexDA = np.load('../data/results/results/distributed_in_time_rh3_200matrixVensembleSplit40state.npz')['msexDA']
# print(msexDA)

# time_arr=[]
# for i in range(10,101,10):
#     time = int(np.load('../data/results/results/matrixVensembleSplit'+str(i)+'elapsed.npy'))
#     time_arr.append(time)
# time_arr.append(int(np.load('../data/results/results/matrixVprec145elapsed.npz')['arr_0']))
#
# print(time_arr)