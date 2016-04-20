```python
>>> import numpy as np
>>> import scipy as sp
>>> import holoviews as hv
>>> hv.notebook_extension()
<IPython.core.display.HTML object>
```

```python
>>> L = 10
>>> nx = 101
>>> psi = np.zeros(nx)
>>> h_bar = 1
>>> mass= 1/2
>>> a = L/(nx-1) # space step
>>> h = 1 # time step
>>> x = np.linspace(0, L, nx)
...
>>> # Potential Well
... V0 = 1
>>> Vwidth = 1
>>> V = np.piecewise(x, [x < L/2 - Vwidth/2, x >= L/2 - Vwidth/2, x > L/2 + Vwidth/2], [0, V0, 0])
...
>>> # Initial wave packet
... E = 1
>>> sigma0 = L/10
>>> p = 0
>>> psi = np.exp(-(x - L/4)**2/(2*sigma0**2)) * np.exp(1j*p*x)
>>> psi /= np.linalg.norm(psi)
...
>>> # Crank Nicolson
... a = - 1j*h/(4*a**2)
>>> b = 0.5 + (1j*h)/(2*a**2) + (1j*h*V)/(4)
>>> H = sp.sparse.diags([a, b, a], [-1, 0, 1], shape=(nx, nx))
>>> A = np.linalg.inv()
>>> A.dot(B)
```

```python
>>> print(A)
[[ 0.5 -8.00000000e-04j  0.0 -2.50000000e+01j  0.0 +0.00000000e+00j ...,
   0.0 +0.00000000e+00j  0.0 +0.00000000e+00j  0.0 +0.00000000e+00j]
 [ 0.0 -2.50000000e+01j  0.5 -8.00000000e-04j  0.0 -2.50000000e+01j ...,
   0.0 +0.00000000e+00j  0.0 +0.00000000e+00j  0.0 +0.00000000e+00j]
 [ 0.0 +0.00000000e+00j  0.0 -2.50000000e+01j  0.5 -8.00000000e-04j ...,
   0.0 +0.00000000e+00j  0.0 +0.00000000e+00j  0.0 +0.00000000e+00j]
 ..., 
 [ 0.0 +0.00000000e+00j  0.0 +0.00000000e+00j  0.0 +0.00000000e+00j ...,
   0.5 -8.00000000e-04j  0.0 -2.50000000e+01j  0.0 +0.00000000e+00j]
 [ 0.0 +0.00000000e+00j  0.0 +0.00000000e+00j  0.0 +0.00000000e+00j ...,
   0.0 -2.50000000e+01j  0.5 -8.00000000e-04j  0.0 -2.50000000e+01j]
 [ 0.0 +0.00000000e+00j  0.0 +0.00000000e+00j  0.0 +0.00000000e+00j ...,
   0.0 +0.00000000e+00j  0.0 -2.50000000e+01j  0.5 -8.00000000e-04j]]
```

```python

```

```python

```

```python

```

```python

```

```python

```

```python

```
