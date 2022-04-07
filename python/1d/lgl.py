from scipy import *
import numpy as np

def lglnodes(n,eps=1.0e-14):
    '''
    Python translation of lglnodes.m

    Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde 
    matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
    integration and spectral methods. 

    Parameters
    ----------
    n : integer, requesting an nth-order Gauss-quadrature rule on [-1, 1]

    Returns
    -------
    (nodes, weights) : tuple, representing the quadrature nodes and weights.
                       Note: (n+1) nodes and weights are returned.
            

    Example
    -------
    >>> from lglnodes import *
    >>> (nodes, weights) = lglnodes(3)
    >>> print(str(nodes) + "   " + str(weights))
    [-1.        -0.4472136  0.4472136  1.       ]   [0.16666667 0.83333333 0.83333333 0.16666667]

    Notes
    -----

    Reference on LGL nodes and weights:  
      C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
      in Fluid Dynamics," Section 2.3. Springer-Verlag 1987

    Written by Greg von Winckel - 04/17/2004
        Contact: gregvw@chtm.unm.edu

    Translated and modified into Python by Jacob Schroder - 9/15/2018 
    '''

    

    w = zeros((n+1,))
    x = zeros((n+1,))
    xold = zeros((n+1,))

    # The Legendre Vandermonde Matrix
    P = zeros((n+1,n+1))

    epss = eps

    # Use the Chebyshev-Gauss-Lobatto nodes as the first guess
    for i in range(n+1): 
        x[i] = -cos(pi*i / n)
  
  
    # Compute P using the recursion relation
    # Compute its first and second derivatives and 
    # update x using the Newton-Raphson method.
    
    xold = 2.0
    
    for i in range(100):
        xold = x
       
        P[:,0] = 1.0 
        P[:,1] = x
       
        for k in range(2,n+1):
            P[:,k] = ( (2*k-1)*x*P[:,k-1] - (k-1)*P[:,k-2] ) / k
       
        x = xold - ( x*P[:,n] - P[:,n-1] )/( (n+1)*P[:,n]) 
        
        if (max(abs(x - xold).flatten()) < epss ):
            break 
    
    w = 2.0 / ( (n*(n+1))*(P[:,n]**2))
    
    return x, w
 

#------------------------------------------------------------------------------
# See section on "Gauss-Lobattto rules" here
#    https://en.wikipedia.org/wiki/Gaussian_quadrature
def gauss_lobatto_points(n):
    if n == 2:
        x = [-1,1]
    elif n == 3:
        x = [-1,0,1]
    elif n == 4:
        x = [-1, -1/np.sqrt(5), 1/np.sqrt(5), 1]
    elif n == 5:
        x = [-1, -np.sqrt(3.0/7.0), 0, np.sqrt(3.0/7.0), 1]
    elif n == 6:
        x = [-1, -np.sqrt(1.0/3.0 + 2.0*np.sqrt(7)/21), 
                 -np.sqrt(1.0/3.0 - 2.0*np.sqrt(7)/21), 
                  np.sqrt(1.0/3.0 - 2.0*np.sqrt(7)/21), 
                  np.sqrt(1.0/3.0 + 2.0*np.sqrt(7)/21), 
              1]
    elif n == 7:
        x = [-1, -np.sqrt(5.0/11.0 + (2.0/11)*np.sqrt(5.0/3.0)),
                 -np.sqrt(5.0/11.0 - (2.0/11)*np.sqrt(5.0/3.0)),
                  0,
                  np.sqrt(5.0/11.0 - (2.0/11)*np.sqrt(5.0/3.0)),
                  np.sqrt(5.0/11.0 + (2.0/11)*np.sqrt(5.0/3.0)),
              1]
    else:
        x, w = lglnodes(n-1)
    return np.array(x)
