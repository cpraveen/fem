import numpy as npy


#-----------------------------------------------------------
# writes an geo geometry file
def writeGEO(filename, X, Y, nWK):
    f = open(filename, 'w')
    print 'Writing ', filename
    
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
    f.write('//Joukowski geometry\n')

    #----------#
    # Vertices #
    #----------#
    floatformat = "{:3.16e}"
    
    for i in range(nnode):
        f.write('Point(' + str(i+1) + ') = {' + floatformat.format(Vx[i]) + ', ' + floatformat.format(Vy[i]) + ', 0.0};\n')
      
    #-------#
    # Faces #
    #-------#
    f.write('BSpline(1) = {')
    for i in range(nAf):
        f.write(str(i+1) + ', ')
    f.write(' 1 };\n') #Close the loop

    f.write('BSpline(2) = {')
    for i in range(nAf, nIn-1):
        f.write(str(i+1) + ', ')
    f.write(str(nIn) + ' };\n')

    f.write('BSpline(3) = {')
    for i in range(nIn-1, nOut-1):
        f.write(str(i+1) + ', ')
    f.write(str(i+1) + ', ' + str(nOut) + ', ' + str(nAf+1) + ' };\n') 

    f.close()
