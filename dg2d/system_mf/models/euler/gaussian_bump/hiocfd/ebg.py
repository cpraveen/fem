import numpy as npy


#-----------------------------------------------------------
# writes an ebg geometry file
def writeEBG(filename, X, Y, nWK):
    f = open(filename, 'w')
    print 'Writing ', filename
    
    f.write('nEmbeddedBoundaryFaceGroups\n')
    f.write('3\n')

    # Airfoil
    Vx = X[nWK+1:-nWK,0]
    Vy = Y[nWK+1:-nWK,0]
    nAf = len(Vx)

    # Inflow
    Vx = npy.append(Vx, X[:,-1])
    Vy = npy.append(Vy, Y[:,-1])
    nIn = len(Vx)

    # Outflow Upper
    Vx = npy.append(Vx, X[-1,-2:0:-1])
    Vy = npy.append(Vy, Y[-1,-2:0:-1])

    # Outflow Lower
    Vx = npy.append(Vx, X[0,:-1] )
    Vy = npy.append(Vy, Y[0,:-1])
    nOut = len(Vx)

    nnode = len(Vx)
    f.write('Vertices\n')
    f.write(str(nnode) + '\n')

    #----------#
    # Vertices #
    #----------#
    floatformat = "{:3.16e}"
    
    for i in range(nnode):
        f.write(floatformat.format(Vx[i]) + ' ' + floatformat.format(Vy[i]) + '\n')
      
    #-------#
    # Faces #
    #-------#
    f.write('Faces\n')
    f.write(str(nAf) + '\n')
    for i in range(nAf-1):
        f.write(str(i+1) + ' ' + str(i+2) + '\n')
    f.write(str(nAf) + ' 1\n') #Close the loop

    f.write(str(nIn-nAf-1) + '\n')
    for i in range(nAf, nIn-1):
        f.write(str(i+1) + ' ' + str(i+2) + '\n')

    f.write(str(nOut-nIn+1) + '\n')
    for i in range(nIn-1, nOut-1):
        f.write(str(i+1) + ' ' + str(i+2) + '\n')
    f.write(str(nOut) + ' ' + str(nAf+1) + '\n') 

    f.close()
