xmin, xmax = -1.0, 1.0
def initial_condition(x):
    from numpy import exp,abs,sqrt,log,empty_like
    f = empty_like(x)
    for i,y in enumerate(x):
        if y > -0.8 and y < -0.6:
            f[i] = exp(-log(2.0)*(y+0.7)**2/0.0009);
        elif y > -0.4 and y < -0.2:
            f[i] = 1.0;
        elif y > 0.0 and y < 0.2:
            f[i] = 1.0 - abs(10.0*(y-0.1));
        elif y > 0.4 and y < 0.6:
            f[i] = sqrt(1.0 - 100.0*(y-0.5)**2);
        else:
            f[i] = 0.0;
    return f
