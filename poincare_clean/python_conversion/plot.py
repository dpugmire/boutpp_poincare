import scipy.io
import matplotlib.pyplot as plt


matFiles = []
matFiles.append('./mat_pp/x1y55z1_v3lc-00-250p.mat')
matFiles.append('./mat_pp/x2y55z1_v3lc-00-250p.mat')
matFiles.append('./mat_pp/x3y55z1_v3lc-00-250p.mat')

matFiles.append('./mat_traj/x1y55z1_v3lc-01-250p.mat')
matFiles.append('./mat_traj/x2y55z1_v3lc-01-250p.mat')
matFiles.append('./mat_traj/x3y55z1_v3lc-01-250p.mat')


for f in matFiles :
    print('FILE: ', f)
    data = scipy.io.loadmat(f)
    dataKeys = data.keys()
    
    for k in dataKeys :
        if 'header' in k or 'version' in k or 'globals' in k : continue
        x = data[k]
        print ('  key: ', k, x.shape)


f = './mat_pp/x3y55z1_v3lc-00-250p.mat'
data = scipy.io.loadmat(f)
X = data['v4']
Y = data['v5']
Z = data['v3']

plt.plot(X, Y, label='Data')
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.title('Plot from .mat file')
plt.legend()
plt.grid(True)
plt.show()



matFileX = './mat_pp/x1y55z1_v3lc-00-250p.mat'
matFileY = './mat_pp/x2y55z1_v3lc-00-250p.mat'
matFileY = './mat_pp/x3y55z1_v3lc-00-250p.mat'

dataX = scipy.io.loadmat(matFileX)
dataY = scipy.io.loadmat(matFileY)

traj = scipy.io.loadmat('./mat_traj/x3y55z1_v3lc-01-250p.mat')
print('8888  ', traj.keys())
print('9999  ', traj['v'].shape)


# Inspect the contents of the .mat file
print("Keys in the .mat file:", dataX.keys())
print("Keys in the .mat file:", dataY.keys())

for k in ['v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v8'] :
    x = dataX[k]
    y = dataY[k]
    print('shapes:', x.shape, y.shape)

y = dataY['v2']
print('---- shape: ', y.shape)
valy = y[0]
print(y)
#print(valy[1])


x = []
n = y.shape[0]
for i in range(n) :
    x.append(i)

plt.plot(x, y, label='Data')
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.title('Plot from .mat file')
plt.legend()
plt.grid(True)
plt.show()
