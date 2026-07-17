from __future__ import division
import numpy as npy
import pylab as pyl

#-----------------------------------------------------------
# writes a 2D ascii plot3d grid file
def writePlot2D(filename, X, Y):
    f = open(filename, 'w')
    print('Writing ', filename)
    ni, nj = X.shape
    
    npy.set_printoptions( precision=16, threshold = ni )
    
    f.write('1'+'\n')
    f.write(str(ni) + ' ' + str(nj) + ' 1'+'\n')
    for j in range(nj):
        f.write(npy.array_str(X[:,j])[1:-1]+'\n')
    for j in range(nj):
        f.write(npy.array_str(Y[:,j])[1:-1]+'\n')
    for j in range(nj):
        f.write('0 '*ni+'\n')
               
    f.close()

#-----------------------------------------------------------
# writes a 2D ascii Laballiur grid file
def writeLaballiur(filename, X, Y, nWK):
    f = open(filename, 'w')
    print('Writing ', filename)
    ni, nj = X.shape
    
    npy.set_printoptions( precision=16, threshold = ni )
    
    f.write(str(ni) + ' ' + str(nj) + ' ' + str(nWK) + ' ' + str(ni-nWK+1) + '\n')
    for j in range(nj):
        f.write(npy.array_str(X[:,j])[1:-1]+'\n')
    for j in range(nj):
        f.write(npy.array_str(Y[:,j])[1:-1]+'\n')


    f.close()

#-----------------------------------------------------------
def Write3DArray(f,X):	
    ni, nj = X.shape
    
    for j in range(nj):
        f.write(npy.array_str(X[:,j])[1:-1]+'\n')

#-----------------------------------------------------------
# writes a 3D ascii plot3d grid file
def writePlot3D(filename, X, Y):
    f = open(filename, 'w')
    print('Writing ', filename)
    ni, nj = X.shape; nk = 2
    
    npy.set_printoptions( precision=16, threshold = ni )
    
    f.write('1'+'\n')
    f.write(str(ni) + ' ' + str(nj) + ' ' + str(nk) + '\n')
    Write3DArray(f,X)
    Write3DArray(f,X)
    Write3DArray(f,Y)
    Write3DArray(f,Y)
    Z = npy.zeros(X.shape)
    Write3DArray(f,Z)
    Z = npy.ones(X.shape)
    Write3DArray(f,Z)
               
    f.close()

#-----------------------------------------------------------
# writes a 3D ascii plot3d grid file
def writeOVERFLOW(filename, X, Y):
    f = open(filename, 'w')
    print('Writing ', filename)
   
    # Overflow requires 3 spanwise ndoes

    ni, nj = X.shape; nk = 3
    
    npy.set_printoptions( precision=16, threshold = ni )
    
    f.write('1'+'\n')
    f.write(str(ni) + ' ' + str(nj) + ' ' + str(nk) + '\n')
    Write3DArray(f,X)
    Write3DArray(f,X)
    Write3DArray(f,X)

    Z = npy.ones(X.shape)*0
    Write3DArray(f,Z)
    Z = npy.ones(X.shape)*-0.5
    Write3DArray(f,Z)
    Z = npy.ones(X.shape)*-1
    Write3DArray(f,Z)

    Write3DArray(f,Y)
    Write3DArray(f,Y)
    Write3DArray(f,Y)
               
    f.close()

#-----------------------------------------------------------
# writes a 3D ascii plot3d grid file
def writePlot3Dxz(filename, X, Y):
    f = open(filename, 'w')
    print('Writing ', filename)
   
    # Overflow requires 3 spanwise ndoes

    ni = 2; nj, nk = X.shape; 
    
    npy.set_printoptions( precision=16, threshold = ni*nj )
    
    f.write('1'+'\n')
    f.write(str(ni) + ' ' + str(nj) + ' ' + str(nk) + '\n')
    
    xx = npy.zeros((ni,nj,nk))
    xx[::2,:,:] = X
    xx[1::2,:,:] = X
    
    for k in range(nk):
        #Write3DArray(f,xx[:,:,k])
        f.write(npy.array_str(xx[:,:,k].flatten('F'))[1:-1]+'\n')
    del xx
    
    yy = npy.ones((ni,nj,nk))
    yy[1,:,:] = 0

    for k in range(nk):
        #Write3DArray(f,yy[:,:,k])
        f.write(npy.array_str(yy[:,:,k].flatten('F'))[1:-1]+'\n')
    del yy
    
    zz = npy.zeros((ni,nj,nk))
    zz[::2,:,:] = Y
    zz[1::2,:,:] = Y

    for k in range(nk):
        #Write3DArray(f,zz[:,:,k])
        f.write(npy.array_str(zz[:,:,k].flatten('F'))[1:-1]+'\n')
    del zz
    
    f.close()

#-----------------------------------------------------------
# reads a 2D ascii plot3d grid file assuming an X-Z plane (i.e. OVERFLOW)
def readPlot2D(filename):
    f = open(filename, 'r')
    print('Reading ', filename)
    
    data = f.read()
    f.close()
    data = data.strip()
    data = npy.array(data.split()).astype(npy.float)
    ngrid = int(data[0])
    assert(ngrid == 1)
    ni = int(data[1]); nj = int(data[2]); nk = int(data[3])
    assert(nj == 1)
        
    X = data[4+0*ni*nk:4+1*ni*nk].reshape( (nk,ni) ).T
    Y = data[4+2*ni*nk:4+3*ni*nk].reshape( (nk,ni) ).T
    
    return X, Y
