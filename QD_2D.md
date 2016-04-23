```python
>>> import numpy as np
>>> import scipy as sp
>>> import holoviews as hv
>>> hv.notebook_extension()
<IPython.core.display.HTML object>
```

```python
>>> L = 100
>>> nx = 201
>>> nt = 300
>>> psi = np.zeros(nx)
>>> h_bar = 1
>>> mass= 1/2
>>> a = L/(nx - 1) # space step
>>> h = 1 # time step
>>> x = np.linspace(0, L, nx)
...
>>> # Potential Well
... V0 = -1.5e1
>>> Vwidth = 10 #0.25*L
>>> Vcenter = 0.25*L
>>> V = np.piecewise(x, [x < Vcenter - Vwidth/2, x >= Vcenter - Vwidth/2, x > Vcenter + Vwidth/2], [0, V0, 0])
>>> #V = np.zeros(nx)
...
... # Initial wave packet
... psi = np.zeros((nx, nt), dtype='cfloat')
>>> sigma0 = L/50
>>> k = 1
>>> psi[:, 0] = np.exp(-(x - L/4)**2/(2*sigma0**2)) * np.exp(1j*k*x)
>>> psi[:, 0] /= np.linalg.norm(psi[:, 0])
```

```python
>>> print(psi.shape)
>>> (hv.Curve(np.absolute(psi[:, 0]), label='Wave function') *
>>>  hv.Curve(V, label='Potential'))
(201, 300)
:Overlay
   .Curve.Wave_function :Curve   [x]   (y)
   .Curve.Potential     :Curve   [x]   (y)
```

```python
>>> # Crank Nicolson
... b = h_bar*h/(4*mass*a**2)
>>> H = sp.sparse.diags([-b, h*V/(2*np.pi) + 2*b, -b], [-1, 0, 1], shape=(nx, nx), dtype='cfloat').todense()
>>> print(np.isfinite(H))
>>> A = sp.linalg.inv(np.eye(nx) + 1j*H)
>>> B = np.eye(nx) - 1j*H
>>> C = np.dot(B, A)
[[ True  True  True ...,  True  True  True]
 [ True  True  True ...,  True  True  True]
 [ True  True  True ...,  True  True  True]
 ..., 
 [ True  True  True ...,  True  True  True]
 [ True  True  True ...,  True  True  True]
 [ True  True  True ...,  True  True  True]]
```

```python
>>> for i in range(nt-1):
...     psi[:, i+1] = np.dot(C, psi[:, i])
```

```python
>>> print(psi.max())
(0.436627762768+0.0402625945862j)
```

```python

```

```python
>>> %output holomap='scrubber'
>>> %output max_frames=100000
>>> hv.HoloMap([(
...             i*h, hv.Curve(np.abs(psi[:, i])**2)) # * hv.Curve(V))
...             for i in range(nt)], kdims = ["Time"])# * hv.Curve(V)
>>> #hv.HoloMap(np.absolute(psi), kdims=['Time'], label='Wave function')
b':HoloMap   [Time]\n   :Curve   [x]   (y)'
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
