import numpy as np
import matplotlib.pyplot as plt

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
        print('Not implemented')
        exit()
    return np.array(x)

# Plot the points
for n in range(2,8):
    plt.plot(gauss_lobatto_points(n), (10-n)*np.ones(n), 'o-')
plt.yticks([])
plt.title('Gauss-Lobatto points')
plt.show()
