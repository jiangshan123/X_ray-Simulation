import numpy as np
from scipy.optimize import leastsq
import pylab as pl

def func(x, p):
    """
    A*sin(2*pi*k*x + theta)
    """
    A, k, theta = p
    return A*np.sin(2*np.pi*k*x+theta)

def residuals(p, y, x):
    """
    the difference between real data and fitting data, p is the fitting index
    """
    return y - func(x, p)

x = np.linspace(0, -2*np.pi, 100)
A, k, theta = 10, 0.34, np.pi/6 #
y0 = func(x, [A, k, theta]) # real data
y1 = y0 + 2 * np.random.randn(len(x)) # adding noise

p0 = [7, 0.2, 0] # first fitting

# Using leastsq to fit
# residuals to calculate error
# p0 is the initial value
# args is real data
plsq = leastsq(residuals, p0, args=(y1, x))

print (u"Real parameters:", [A, k, theta])
print (u"Fitting parameters", plsq[0]) # 实验数据拟合后的参数

pl.plot(x, y0, label=u"Real data")
pl.plot(x, y1, label=u"Data with noise")
pl.plot(x, func(x, plsq[0]), label=u"Fitting data")
pl.legend()
pl.show()