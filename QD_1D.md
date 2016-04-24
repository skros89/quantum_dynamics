## 1D Quantum Dynamics

```python
>>> import numpy as np
>>> import scipy as sp
>>> import holoviews as hv
>>> hv.notebook_extension()
<IPython.core.display.HTML object>
```

```python
>>> L = 12.5 # x-domain from 0 to L
>>> nx = 2801
>>> nt = 3000
>>> a = L/(nx - 1) # space step
>>> h = 2.5e-5 # time step
>>> x = np.linspace(0, L, nx)
>>> dk = 2*np.pi/L
>>> k0 = 30
>>> k = np.fft.fftfreq(nx, a)
...
>>> # Potential Well
... V0 = 1500
>>> Vwidth = 0.01*L
>>> Vcenter = 0.6*L
>>> V = np.piecewise(x, [x < Vcenter - Vwidth/2, x >= Vcenter - Vwidth/2, x > Vcenter + Vwidth/2], [0, V0, 0])
>>> #V = np.zeros(nx)
...
... # Initial wave packet
... psi = np.zeros((nx, nt), dtype='cfloat')
>>> sigma0 = 0.1
>>> psi[:, 0] = np.exp(-(x - L/2)**2/(2*sigma0**2)) * np.exp(1j*k0*x)
>>> psi[:, 0] /= np.linalg.norm(psi[:, 0])
>>> %opts Curve norm{+axiswise}
>>> (hv.Curve(np.absolute(psi[:, 0]), label='Wave function') + hv.Curve(V, label='Potential'))
:Layout
   .Curve.Wave_function :Curve   [x]   (y)
   .Curve.Potential     :Curve   [x]   (y)
```

## Crank Nicolson

```python
>>> # Crank Nicolson
... b = h/(2*a**2)
>>> H = sp.sparse.diags([-b, h*V/(2*np.pi) + 2*b, -b], [-1, 0, 1], shape=(nx, nx), dtype='cfloat').todense()
>>> A = sp.linalg.inv(np.eye(nx) + 1j*H)
>>> B = np.eye(nx) - 1j*H
>>> C = np.dot(B, A)
...
>>> for i in range(nt-1):
...     psi[:, i+1] = np.dot(C, psi[:, i])
...
>>> %output holomap='scrubber'
>>> %output max_frames=100000
>>> hv.HoloMap([(i*h, hv.Curve(np.abs(psi[:, i])**2))for i in range(nt)], kdims = ["Time"])
b':HoloMap   [Time]\n   :Curve   [x]   (y)'
```

## Split-operator scheme

```python
>>> nt = 10000
>>> psi = np.zeros((nx, nt), dtype='cfloat')
>>> sigma0 = 0.1
>>> psi[:, 0] = np.exp(-(x - L/2)**2/(2*sigma0**2)) * np.exp(1j*k0*x)
>>> psi[:, 0] /= np.linalg.norm(psi[:, 0])
>>> for i in range(nt-1):
...     psi_p = np.fft.fft(np.exp(-1j*h*V)*psi[:, i])
...     psi[:, i+1] = np.fft.ifft(np.exp(-1j*h*k**2)*psi_p)
...
>>> hv.HoloMap([(i*h, hv.Curve(np.abs(psi[:, i])**2))for i in range(nt)], kdims = ["Time"])
```

```python

```
