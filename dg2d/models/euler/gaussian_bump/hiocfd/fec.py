from __future__ import division

def writeFEC(filename_base, ref, Q, E, V, nLE, NC, nWK, nWB, nr):
    #=========================#
    # Write out the grid file #
    # UBC curved FE format    #
    #                         #
    # Also writes the boundary#
    # shape in GRUMMP .bdry   #
    # format                  #
    #=========================#

    if (Q != 3):
        print("Error: can't write cubic FE data with Q = " + str(Q) + '\n');
        return;
    
    filename = filename_base + '_ref'+str(ref)+'.fec'
    print 'Writing ', filename
    f = open(filename, 'w')

    nelem = E.shape[0];
    nnode = V.shape[0];

    f.write(str(nnode) + ' ' + str(nelem) + ' 4\n');
    
    #----------#
    # Vertices #
    #----------#
    floatformat = "{:3.16e}"
    for i in range(nnode):
        f.write(floatformat.format(V[i,0]) + ' ' +
                floatformat.format(V[i,1]) + '\n')
      
    #----------#
    # Elements #
    #----------#
    # Write as cubic quads
    
    for e in range(nelem):
        f.write('12   ');
        # First write the four corners
        f.write(str(E[e,0]-1)+' ');
        f.write(str(E[e,3]-1)+' ');
        f.write(str(E[e,15]-1)+' ');
        f.write(str(E[e,12]-1)+' ');
        # Now write the edge nodes, CCW from vert 0
        f.write(str(E[e,1]-1)+' ');
        f.write(str(E[e,2]-1)+' ');
        f.write(str(E[e,7]-1)+' ');
        f.write(str(E[e,11]-1)+' ');
        f.write(str(E[e,14]-1)+' ');
        f.write(str(E[e,13]-1)+' ');
        f.write(str(E[e,8]-1)+' ');
        f.write(str(E[e,4]-1)+' ');
        f.write('\n');

    f.close()

    f = open(filename_base + '_ref'+str(ref)+'.bdry', 'w')
    #----------------#
    # Boundary faces #
    #----------------#

    # Number of points, number of curves
    f.write(str(nLE + nWB + 2*nr - 4) + ' 3\n');
    # Airfoil coords
    f.write('# airfoil coordinates ' + str(nLE-1) + ' points\n');
    for i in range(0,nLE-1):
        index = NC[nWK-2+i,0]
        f.write(floatformat.format(V[index,0]) + ' ' +
                floatformat.format(V[index,1]) + '\n')
#        f.write(str(i) + ' ' + floatformat.format(V[index,0]) + ' ' +
#                floatformat.format(V[index,1]) + '\n')

      
    # Farfield inflow
    f.write('#inflow ' + str(nWB) + ' points\n')
    for i in range(nWB):
        index = NC[i, nr-1] - 1
        f.write(floatformat.format(V[index,0]) + ' ' +
                floatformat.format(V[index,1]) + '\n')
#        f.write(str(i  +nLE-1) + ' ' + floatformat.format(V[index,0]) + ' ' +
#                floatformat.format(V[index,1]) + '\n')

    # Farfield Outflow
    nb = int(2*(nr-1));
    f.write('#outflow ' + str(nb) + ' points \n')
    for i in range(1,nr-1):
        index = NC[0,nr - i - 1] - 1
        f.write(floatformat.format(V[index,0]) + ' ' +
                floatformat.format(V[index,1]) + '\n')
#        f.write(str(i + nLE-1 + nWB-1) + ' ' + floatformat.format(V[index,0]) + ' ' +
#                floatformat.format(V[index,1]) + '\n')
    f.write('#outflow part 2\n')
    for i in range(nr-1):
        index = NC[nWB-1,i] - 1
        f.write(floatformat.format(V[index,0]) + ' ' +
                floatformat.format(V[index,1]) + '\n')
#        f.write(str(i + nLE-1 + nWB-1 + nr) + ' ' + floatformat.format(V[index,0]) + ' ' +
#                floatformat.format(V[index,1]) + '\n')

    # Airfoil surface
    f.write('spline r 1 b 1 ' + str(nLE) + ' ');
    for i in range(0, nLE-1):
        f.write(str(i) + ' ')
    f.write(' 0\n')
    
    # Inflow part of domain
    f.write('spline b 5 r 1 ' + str(nWB) + ' ')
    for i in range(nWB):
        f.write(str(i+nLE-1) + ' ')
    f.write('\n')
      
    # Outflow end of the domain
    f.write('spline r 1 b 5 ' + str(nb+1) + ' ' + str(nLE-1) + ' ')
    for i in range(nb-1):
        f.write(str(i + nLE + nWB - 1) + ' ')
    f.write(str(nLE-1 + nWB-1) + '\n')
    f.close()

    return
