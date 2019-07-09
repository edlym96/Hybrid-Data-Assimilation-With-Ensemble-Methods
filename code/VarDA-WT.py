import numpy as np
from scipy.optimize import minimize
import random
import matplotlib.pyplot as plt
import time

from numpy.linalg import inv
from numpy import linalg as LA


import math
from scipy.sparse.linalg import svds


import sys
sys.path.append('fluidity-master')
import vtktools


ntime = 5000
ntimelocal = 500
#ntime = 500
#ntime = 30






'''


uvwTot = np.array([])
#for i in range(100,ntimelocal):
for i in range(ntimelocal):

        filename = '/media/rarcucci/Maxtor/For-Rosella-Dec18/Projected_Normal_'+str(i+1)+'.vtu'
        ug=vtktools.vtu(filename)
        ug.GetFieldNames()
        uvwVec = ug.GetScalarField('Tracer1')
        dimuvw=len(uvwVec)
        n=len(uvwVec)
        uvwTot= np.append(uvwTot,uvwVec)


for i in range(ntimelocal):

        filename = '/media/rarcucci/Maxtor/For-Rosella-Dec18/Projected_Normal_'+str(i+1)+'.vtu'
        ug=vtktools.vtu(filename)
        ug.GetFieldNames()
        uvwVec = ug.GetScalarField('Tracer1')
        dimuvw=len(uvwVec)
        n=len(uvwVec)
        uvwTot= np.append(uvwTot,uvwVec)

for i in range(ntimelocal):

        filename = '/media/rarcucci/Maxtor/For-Rosella-Dec18/Projected_Normal_'+str(i+1)+'.vtu'
        ug=vtktools.vtu(filename)
        ug.GetFieldNames()
        uvwVec = ug.GetScalarField('Tracer1')
        dimuvw=len(uvwVec)
        n=len(uvwVec)
        uvwTot= np.append(uvwTot,uvwVec)


for i in range(ntimelocal):

        filename = '/media/rarcucci/Maxtor/For-Rosella-Dec18/Projected_Normal_'+str(i+1)+'.vtu'
        ug=vtktools.vtu(filename)
        ug.GetFieldNames()
        uvwVec = ug.GetScalarField('Tracer1')
        dimuvw=len(uvwVec)
        n=len(uvwVec)
        uvwTot= np.append(uvwTot,uvwVec)

for i in range(ntimelocal):

        filename = '/media/rarcucci/Maxtor/For-Rosella-Dec18/Projected_Normal_'+str(i+1)+'.vtu'
        ug=vtktools.vtu(filename)
        ug.GetFieldNames()
        uvwVec = ug.GetScalarField('Tracer1')
        dimuvw=len(uvwVec)
        n=len(uvwVec)
        uvwTot= np.append(uvwTot,uvwVec)


for i in range(ntimelocal):

        filename = '/media/rarcucci/Maxtor/For-Rosella-Dec18/Projected_Normal_'+str(i+1)+'.vtu'
        ug=vtktools.vtu(filename)
        ug.GetFieldNames()
        uvwVec = ug.GetScalarField('Tracer1')
        dimuvw=len(uvwVec)
        n=len(uvwVec)
        uvwTot= np.append(uvwTot,uvwVec)

for i in range(ntimelocal):

        filename = '/media/rarcucci/Maxtor/For-Rosella-Dec18/Projected_Normal_'+str(i+1)+'.vtu'
        ug=vtktools.vtu(filename)
        ug.GetFieldNames()
        uvwVec = ug.GetScalarField('Tracer1')
        dimuvw=len(uvwVec)
        n=len(uvwVec)
        uvwTot= np.append(uvwTot,uvwVec)


for i in range(ntimelocal):

        filename = '/media/rarcucci/Maxtor/For-Rosella-Dec18/Projected_Normal_'+str(i+1)+'.vtu'
        ug=vtktools.vtu(filename)
        ug.GetFieldNames()
        uvwVec = ug.GetScalarField('Tracer1')
        dimuvw=len(uvwVec)
        n=len(uvwVec)
        uvwTot= np.append(uvwTot,uvwVec)


for i in range(ntimelocal):

        filename = '/media/rarcucci/Maxtor/For-Rosella-Dec18/Projected_Normal_'+str(i+1)+'.vtu'
        ug=vtktools.vtu(filename)
        ug.GetFieldNames()
        uvwVec = ug.GetScalarField('Tracer1')
        dimuvw=len(uvwVec)
        n=len(uvwVec)
        uvwTot= np.append(uvwTot,uvwVec)


for i in range(ntimelocal):

        filename = '/media/rarcucci/Maxtor/For-Rosella-Dec18/Projected_Normal_'+str(i+1)+'.vtu'
        ug=vtktools.vtu(filename)
        ug.GetFieldNames()
        uvwVec = ug.GetScalarField('Tracer1')
        dimuvw=len(uvwVec)
        n=len(uvwVec)
        uvwTot= np.append(uvwTot,uvwVec)




#ntime = 400
#ntimelocal = 400





dimuvwTot = len(uvwTot)

m = np.array([])
m = np.zeros(dimuvwTot/ntime)
for j in range(n):
        for i in range(1,ntime+1):
                m[j] = np.copy(m[j] + uvwTot[j+(i-1)*n])
        m[j] = m[j]/ntime



err = np.array([])
err = np.zeros(dimuvwTot)

for j in range(n):
        for i in range(1,ntime+1):
                err[j+(i-1)*j] = abs(uvwTot[j+(i-1)*j]-m[j])

W =  np.transpose(np.reshape(err, (ntime,n)))

'''













'''

mid = int(np.floor(ntime/2))

U, s, UT = svds(W, k=mid)

print('il primo svd calcolato!')


st = np.sort(s)[::-1]

s1= st[0]

ref=math.sqrt(s1)
print "ref", ref


trnc=2

while st[trnc] >= ref and trnc < mid-1 :
        trnc=trnc+1
        print "trnc", trnc
        print "st[trnc]", st[trnc]

###########################################

'''

'''

trnc=ntime-1

print('value of trnc')
print(trnc)


Utrunc, strunc, Wtrunc = svds(W, k=trnc)
X = Utrunc.dot(np.diag(np.sqrt(strunc)))
np.savetxt("/home/rarcucci/VarDA/WindTunnel/Results/matrixVprecWT1000"+str(trnc)+".txt", X)


V = X.copy()

'''

# EXACTLY THE SAME AS VARDA-3DTRACERS+COVARIANCE BUT DOESNT SAVE RESULTS

trnc=ntime-1

lam = 1e-60

#lam = 1


V = np.loadtxt("/home/rarcucci/WindTunnel/Results/matrixVprecWT1000"+str(trnc)+".txt", usecols=range(trnc))


maxV = V.max()
minV = V.min()

print 'maxV' , maxV , '\n'
print 'minV' , minV , '\n'








#elapsedvec = np.array([])
#errxBvec = np.array([])
#errxDAvec = np.array([])

#for j in range(300,ntime-1):
ugg=vtktools.vtu('/home/rarcucci/dataWT/TracerWT400.vtu')
ugg.GetFieldNames()
uvwVecobs = ugg.GetScalarField('n2_TracerFluidity_WT')
#uvwVecobs = 100000 * uvwVecobs
nstobs = len(uvwVecobs)
ug=vtktools.vtu('/home/rarcucci/dataWT/For-Rosella-Dec18/Projected_Normal_400.vtu')
ug.GetFieldNames()
uvwVec = ug.GetScalarField('Tracer1')
#uvwVec = 100000 * uvwVec
nst = len(uvwVec)
pos=ug.GetLocations()
z=pos[:,2]
n = len(uvwVec)


diff = np.absolute(uvwVec - uvwVecobs)



m = trnc
xB = uvwVec.copy()
y = uvwVecobs.copy()
R = lam * 0.9
#R = lam


#x0 = 0.0000001 * np.ones(n)
x0 = uvwVec
#V = np.absolute(V)

maxV = V.max()
minV = V.min()

print 'maxV-abs' , maxV , '\n'
print 'minV-abs' , minV , '\n'


Vin = np.linalg.pinv(V)
v0 = np.dot(Vin,x0)
VT = np.transpose(V)
HxB = xB.copy()
d = np.subtract(y,HxB)

# Cost function J
def J(v):
        vT = np.transpose(v)
        vTv = np.dot(vT,v)
        Vv = np.dot(V,v)
        Jmis = np.subtract(Vv,d)
#       invR = 1/R
        invR = 1e+60
        JmisT = np.transpose(Jmis)
        RJmis = JmisT.copy()
        J1 = invR*np.dot(Jmis,RJmis)
        Jv = (vTv + J1) / 2
        return Jv

# Gradient of J
def gradJ(v):
        Vv = np.dot(V,v)
        Jmis = np.subtract(Vv,d)
#       invR = 1/R
        invR = 1e+60
        g1 = Jmis.copy()
        VT = np.transpose(V)
        g2 = np.dot(VT,g1)
        gg2 = np.multiply(invR , g2)
        ggJ = v + gg2
        return ggJ

# Compute the minimum
t = time.time()

res = minimize(J, v0, method='L-BFGS-B', jac=gradJ,
                options={'disp': True})


vDA = np.array([])
vDA = res.x
deltaxDA = np.dot(V,vDA)
xDA = xB + deltaxDA

elapsed = time.time() - t
#print 'elapsed' , elapsed , '\n'

errxB = y - xB
MSExb = LA.norm(errxB, 2)/LA.norm(y, 2)
print 'L2 norm of the background error' , MSExb , '\n'

errxDA = y - xDA
MSExDA = LA.norm(errxDA, 2)/LA.norm(y, 2)
print 'L2 norm of the error in DA solution' , MSExDA , '\n'


#elapsedvec = np.append(elapsedvec,elapsed)
#errxBvec = np.append(errxBvec,MSExb)
#errxDAvec = np.append(errxDAvec,MSExDA)



#np.savetxt("/home/rarcucci/VarDA/TestxDA/Results/elapsedvecTime500-tracers-George.txt",elapsedvec)
#np.savetxt("/home/rarcucci/VarDA/TestxDA/Results/errxBvecTime500-tracers-George.txt",errxBvec)
#np.savetxt("/home/rarcucci/VarDA/TestxDA/Results/errxDAvecTime500-tracers-George.txt",errxDAvec)







abserrxB = np.absolute(errxB)

#vabserrxB = np.reshape(abserrxB, (nst,1))


ug.AddScalarField('u_0^M - u_C', abserrxB)
ug.Write('/home/rarcucci/WindTunnel/Results/abserruM400-new.vtu')

abserrxDA = np.absolute(errxDA)

#vabserrxDA = np.reshape(abserrxDA, (nst,3))

ug.AddScalarField('u^DA - u_C', abserrxDA)
ug.Write('/home/rarcucci/WindTunnel/Results/abserruDA400-new.vtu')


#deltaxDAM = np.reshape(deltaxDA, (nst,3))

#ug.AddVectorField('deltaxDA', deltaxDAM)
#ug.Write('/home/rarcucci/VarDA/TestxDA/Results/deltaxDAPOD-3DVel.vtu')

#xDAM = np.reshape(xDA, (nst,3))

ug.AddScalarField('uDA', xDA)
ug.Write('/home/rarcucci/WindTunnel/Results/uDA400-new.vtu')



#xBM = np.reshape(xB, (nst,3))

ug.AddScalarField('uM', xB)
ug.Write('/home/rarcucci/WindTunnel/Results/uM400-new.vtu')



#yM= np.reshape(y, (nst,3))

ug.AddScalarField('y', y)
ug.Write('/home/rarcucci/WindTunnel/Results/y400-new.vtu')


ug.AddScalarField('Misfit d', d)
ug.Write('/home/rarcucci/WindTunnel/Results/d400-new.vtu')

ug.AddScalarField('Misfit diff', diff)
ug.Write('/home/rarcucci/WindTunnel/Results/diff400-new.vtu')



'''
ug=vtktools.vtu('/du-4tb/Rossella/24POD/LSBU_POD_4.vtu')
ug.GetFieldNames()
#p=ug.GetScalarField('Pressure')
uvwPOD = ug.GetVectorField('Velocity')
nstPOD = len(uvwPOD)
uvwVecPOD = np.ravel(uvwPOD, order='F')

xDAPOD = uvwVecPOD + deltaxDA

xDAPODM = np.reshape(xDAPOD, (nst,3))

ug.AddVectorField('xDAPOD', xDAPODM)
ug.Write('/home/rarcucci/VarDA/TestxDA/Results/xDAPOD4-3DVel.vtu')

'''
