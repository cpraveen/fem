"""
Solve scalar conservation law with periodic bc
To get help, type
    python dg.py -h
"""
import argparse
from dgsolver import *

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-pde', choices=('linear','varadv','burger'), help='PDE', 
                    default='linear')
parser.add_argument('-ncell', type=int, help='No of cells', default=50)
parser.add_argument('-nrefine', type=int, help='No of refinements', default=0)
parser.add_argument('-degree', type=int, help='Polynomial degree', default=1)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-plot_freq', type=int, help='Frequency to plot solution', 
                    default=1)
parser.add_argument('-ic', choices=('sin2pi','sin4pi','gauss','hat','mult'),
                    help='Initial condition', default='sin2pi')
parser.add_argument('-limit', choices=('no','yes'), help='Apply limiter', 
                    default='no')
parser.add_argument('-tvbM', type=float, help='TVB M parameter', default=0.0)
parser.add_argument('-compute_error', choices=('no','yes'), 
                    help='Compute error norm', default='no')
parser.add_argument('-num_flux', choices=('central','upwind','roe','godunov'),
                    help='Numerical flux', default='upwind')
args = parser.parse_args()

if args.nrefine == 0: # Run on a single grid, plot solution
    dg = DG(args)
    dg.run()
    plt.show()
else:                 # Perform grid refinement, disable plotting
    args.plot_freq = 0
    args.compute_error = 'yes'
    error = np.empty(args.nrefine)
    for n in range(args.nrefine):
        dg = DG(args)
        nc, dx, error[n] = dg.run()
        if n == 0:
            print('ncell, dx, error, rate = %5d %12.4e %12.4e' % 
                  (nc, dx, error[n]))
        else:
            rate = np.log(error[n-1]/error[n]) / np.log(2)
            print('ncell, dx, error, rate = %5d %12.4e %12.4e %8.3f' % 
                  (nc, dx, error[n], rate))
        args.ncell *= 2
