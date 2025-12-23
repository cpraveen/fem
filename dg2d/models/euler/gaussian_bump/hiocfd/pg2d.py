from __future__ import division

def writePG2D(filename_base, ref, Q, TriFlag, E, V, nLE, NC, nWK, nWB, nr):
    #=========================#
    # Write out the grid file #
    #=========================#

    assert Q == 1

    filename = filename_base + ('_tri' if TriFlag else '_quad') + '_ref'+str(ref)+ '_Q'+str(Q)+'.pg2d'
    print 'Writing ', filename
    f = open(filename, 'w')
     
    nelem = E.shape[0];
    neli = int((NC.shape[0]-1)/Q)
    nelj = int((NC.shape[1]-1)/Q)

    nAf = int((nLE-1)/Q)
    nFFi = int((nWB-1)/Q)
    nFFo = int((nr-1)/Q)
    
    fac = 2 if TriFlag else 1
    
    nedges = (neli+1)*(nelj+1)*fac + (nelem if TriFlag else 0) - int((nWK-1)/Q)
    
    nbedges = nAf + nFFi + 2*nFFo
    
    ntri = fac*nelem + nbedges
    
    nnode = V.shape[0];
    f.write(str(nnode) + ' ' + str(nedges) + ' ' + str(ntri) + ' ' + str(nbedges) + '\n') #dim nNodes negrp nbfgrp
     
    #----------------#
    # Boundary edges #
    #----------------#
    
    nbc = 1
    # Airfoil
    for i in range(nAf):
        f.write(str(NC[nWK-1+Q*i,0]) + ' ' + str(NC[nWK-1+Q*(i+1),0]) + ' ' + int(nWK-1)/Q+i*fac+1 + ' ' + ntri+nbc + '\n')
        nbc+=1
    
      
    # Farfield inflow
    for i in range(nFFi):
        f.write(str(NC[Q*i,nr-1]) + ' ' + str(NC[Q*(i+1),nr-1]) + ' ' + (neli-1)*nelj*fac+i*fac+1 + ' ' + ntri+nbc + '\n')
        nbc+=1
        
    # Farfield Outflow
    for i in range(nFFo):
        f.write(str(NC[0,Q*i]) + ' ' + str(NC[0,Q*(i+1)]) + ' ' + ntri+nbc + '\n')
        nbc+=1

    for i in range(nFFo):
        f.write(str(NC[nWB-1,Q*i]) + ' ' + str(NC[nWB-1,Q*(i+1)]) + ' ' + ntri+nbc + '\n')
        nbc+=1

    #----------------#
    # Interior Edges #
    #----------------#
    
    #Wake edges
    j = 0
    for i in range(int((nWK-1)/Q)):
        cell1 = i+1
        cell2 = neli - i+1
        f.write(str(NC[Q*i,Q*j]) + ' ' + str(NC[Q*i,Q*(j+1)]) + ' ' + str(cell1) + ' ' + str(cell2) + '\n')

    # i-constant edges
    for j in range(nelj):
        for i in range(neli-1):
            cell1 = i   + neli*j + 1
            cell2 = i+1 + neli*j + 1
            f.write(str(NC[Q*i,Q*j]) + ' ' + str(NC[Q*i,Q*(j+1)]) + ' ' + str(cell1) + ' ' + str(cell2) + '\n')

    # j-constant edges
    for j in range(nelj-1):
        for i in range(neli):
            cell1 = i + neli*j     + 1
            cell2 = i + neli*(j+1) + 1
            f.write(str(NC[Q*i,Q*j]) + ' ' + str(NC[Q*(i+1),Q*j]) + ' ' + str(cell1) + ' ' + str(cell2) + '\n')


    #----------#
    # Vertices #
    #----------#
    floatformat = "{:3.16e}"
    f.write('ndoes')
    for i in range(nnode):
        f.write(floatformat.format(V[i,0]) + ' ' + floatformat.format(V[i,1]) + '\n')

    
    f.close()
    return
