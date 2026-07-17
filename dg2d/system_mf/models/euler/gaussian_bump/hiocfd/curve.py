import numpy as npy

#-----------------------------------------------------------
# writes a curve geometry files
def writeCurve(X, Y, nWK):
    
    # Airfoil
    Afx = X[nWK:-nWK,0]
    Afy = Y[nWK:-nWK,0]
    nAf = len(Afx)

    # Inflow
    Inx = X[-1::-1,-1]
    Iny = Y[-1::-1,-1]
    nIn = len(Inx)

    # Outflow Upper
    Outx = X[0,-1:0:-1]
    Outy = Y[0,-1:0:-1]

    # Outflow Lower
    Outx = npy.append(Outx, X[-1,:] )
    Outy = npy.append(Outy, Y[-1,:])
    nOut = len(Outx)

    filename = 'joukowski_af.csv'
    Af = open(filename, 'w')
    print 'Writing ', filename

    filename = 'joukowski_in.csv'
    In = open(filename, 'w')
    print 'Writing ', filename

    filename = 'joukowski_out.csv'
    Out = open(filename, 'w')
    print 'Writing ', filename

    #----------#
    # Vertices #
    #----------#
    floatformat = "{:3.16e}"
    
    for i in range(nAf):
        Af.write(floatformat.format(Afx[i]) + ' ' + floatformat.format(Afy[i]) + ' 0' + '\n')
    Af.close()

    for i in range(nIn):
        In.write(floatformat.format(Inx[i]) + ' ' + floatformat.format(Iny[i]) + ' 0' + '\n')
    In.close()

    for i in range(nOut):
        Out.write(floatformat.format(Outx[i]) + ' ' + floatformat.format(Outy[i]) + ' 0' + '\n')
    Out.close()
