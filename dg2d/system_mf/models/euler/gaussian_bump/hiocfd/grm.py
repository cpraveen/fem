from __future__ import division

def writeGRM(filename_base, ref, Q, TriFlag, E, V, nLE, NC, nWK, nWB, nr):
    #=========================#
    # Write out the grid file #
    #=========================#

    filename = filename_base + '_ref'+str(ref)+ '_Q'+str(Q)+'.grm'
    print 'Writing ', filename
    f = open(filename, 'w')
     
    nelem = E.shape[0];
    
    fac = 2 if TriFlag else 1
    
    nnode = V.shape[0];
    f.write('2 ' + str(nnode) + ' 1 3\n') #dim nNodes negrp nbfgrp
     
    #----------#
    # Vertices #
    #----------#
    floatformat = "{:3.16e}"
    
    for i in range(nnode):
        f.write(floatformat.format(V[i,0]) + ' ' + floatformat.format(V[i,1]) + '\n')
      
    #----------------#
    # Boundary faces #
    #----------------#
        
    # Airfoil
    nb = int((nLE-1)/Q)
    f.write(str(nb) + '\n');
    f.write('PXE_Shape_Edge\n')
    for i in range(int((nLE-1)/Q)):
        f.write(str(NC[nWK-1+Q*i,0]) + ' ' + str(NC[nWK-1+Q*(i+1),0]) + '\n')
    
      
    # Farfield inflow
    f.write(str(int((nWB-1)/Q)) + '\n')
    f.write('PXE_Shape_Edge\n')
    for i in range(int((nWB-1)/Q)):
        f.write(str(NC[Q*i,nr-1]) + ' ' + str(NC[Q*(i+1),nr-1]) + '\n')
      
    # Farfield Outflow
    nb = int(2*(nr-1)/Q);
    f.write(str(nb) + '\n')
    f.write('PXE_Shape_Edge\n')
    for i in range(int((nr-1)/Q)):
        f.write(str(NC[0,Q*i]) + ' ' + str(NC[0,Q*(i+1)]) + '\n')

    for i in range(int((nr-1)/Q)):
        f.write(str(NC[nWB-1,Q*i]) + ' ' + str(NC[nWB-1,Q*(i+1)]) + '\n')
    
      
      
    #----------#
    # Elements #
    #----------#
      
    if TriFlag:
        f.write(str(2*nelem)+'\n')
        f.write(str(Q)+'\n')
        f.write('PXE_Shape_Triangle\n')
        f.write('UniformNodeDistribution\n')

        ni = int((NC.shape[0]-1)/Q)
        nj = int((NC.shape[1]-1)/Q)
        
        if Q == 1:
            nodemap = (0,  1, 
                       2, -1)
            nodemap2= (0, 1, 
                      -1, 2)
        if Q == 2:
            nodemap = (0,  5,  1, 
                       4,  3, -1,   
                       2, -1, -1)
            
            nodemap2= ( 0,  5,  1, 
                       -1,  4,  3,   
                       -1, -1,  2)
    
        if Q == 3:
            nodemap = ( 0,  7,  8,  1, 
                        6,  9,  3, -1, 
                        5,  4, -1, -1, 
                        2, -1, -1, -1)
            nodemap2= ( 0,  7,  8,  1, 
                       -1,  6,  9,  3, 
                       -1, -1,  5,  4, 
                       -1, -1, -1,  2)
        if Q == 4:
            nodemap = ( 0,  9, 10, 11,  1, 
                        8, 12, 13,  3, -1, 
                        7, 14,  4, -1, -1, 
                        6,  5, -1, -1, -1, 
                        2, -1, -1, -1, -1, )
            nodemap2= ( 0,  9, 10, 11,  1, 
                       -1,  8, 12, 13,  3, 
                       -1, -1,  7, 14,  4,
                       -1, -1, -1,  6,  5, 
                       -1, -1, -1, -1,  2,)
        
        #Invert the map
        nodemapinv  = []
        for k in range(int((Q+1)*(Q+2)/2)):
            j = 0
            while nodemap[j] != k: j += 1
            nodemapinv.append(j)
    
        nodemapinv2  = []
        for k in range(int((Q+1)*(Q+2)/2)):
            j = 0
            while nodemap2[j] != k: j += 1
            nodemapinv2.append(j)

        for j in range(nj):
            for i in range(int(ni/2)):
                e = i + ni*j
                
                #Write nodes
                for k in range(int((Q+1)*(Q+2)/2)):
                    f.write(str(E[e,nodemapinv[k]])+' ')
                f.write('\n')
        
                
                #Write nodes
                for k in range(int((Q+1)*(Q+2)/2)):
                    f.write(str(E[e,(Q+1)*(Q+1)-1-nodemapinv[k]])+' ')
                f.write('\n')
                
            for i in range(int(ni/2),ni):
                e = i + ni*j
                
                #Write nodes
                for k in range(int((Q+1)*(Q+2)/2)):
                    f.write(str(E[e,nodemapinv2[k]])+' ')
                f.write('\n')
                        
                #Write nodes
                for k in range(int((Q+1)*(Q+2)/2)):
                    f.write(str(E[e,(Q+1)*(Q+1)-1-nodemapinv2[k]])+' ')
                f.write('\n')

    else:
        f.write(str(nelem)+'\n')
        f.write(str(Q)+'\n')
        f.write('PXE_Shape_Quad\n')
        f.write('UniformNodeDistribution\n')
        #Write nodes
        for e in range(nelem):
            for k in range((Q+1)*(Q+1)):
                f.write(str(E[e,k])+' ')
            f.write('\n')
      
    f.close()
    return
