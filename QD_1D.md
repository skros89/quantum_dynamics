```python
>>> import numpy as np
>>> import scipy as sp
>>> import holoviews as hv
>>> hv.notebook_extension()
<IPython.core.display.HTML object>
```

## Crank Nicolson

```python
>>> L = 1
>>> nx = 2000
>>> nt = 400
>>> psi = np.zeros(nx)
>>> h_bar = 1
>>> mass= 1/2
>>> a = L/(nx - 1) # space step
>>> h = 1e-5 # time step
>>> x = np.linspace(0, L, nx)
...
>>> # Potential Well
... V0 = 10000
>>> Vwidth = 0.25*L
>>> Vcenter = 0.5*L
>>> V = np.piecewise(x, [x < Vcenter - Vwidth/2, x >= Vcenter - Vwidth/2, x > Vcenter + Vwidth/2], [0, V0, 0])
>>> V = np.zeros(nx)
...
>>> # Initial wave packet
... psi = np.zeros((nx, nt), dtype='cfloat')
>>> E = 1
>>> sigma0 = L/50
>>> k = 50
>>> psi[:, 0] = np.exp(-(x - L/8)**2/(2*sigma0**2)) * np.exp(1j*k*x)
>>> psi[:, 0] /= np.linalg.norm(psi[:, 0])
...
>>> (hv.Curve(np.absolute(psi[:, 0]), label='Wave function') *
>>>  hv.Curve(V, label='Potential'))
...
>>> # Crank Nicolson
... b = h_bar*h/(4*mass*a**2)
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
>>> hv.HoloMap([(
...             i*h, hv.Curve(np.real(psi[:, i])**2)) # * hv.Curve(V))
...             for i in range(nt)], kdims = ["Time"])
>>> #hv.HoloMap(np.absolute(psi), kdims=['Time'], label='Wave function')
```

## Trotter-Suzuki

```python
>>> L = 1
>>> nx = 2000
>>> nt = 400
>>> psi = np.zeros(nx)
>>> h_bar = 1
>>> mass= 1/2
>>> a = L/(nx - 1) # space step
>>> h = 1e-5 # time step
>>> x = np.linspace(0, L, nx)
>>> k = np.fft.fftfreq(nx, a)
...
>>> # Potential Well
... V0 = -3000
>>> Vwidth = 0.25*L
>>> Vcenter = 0.5*L
>>> V = np.piecewise(x, [x < Vcenter - Vwidth/2, x >= Vcenter - Vwidth/2, x > Vcenter + Vwidth/2], [0, V0, 0])
>>> V = np.zeros(nx)
...
>>> # Initial wave packet
... psi = np.zeros((nx, nt), dtype='cfloat')
>>> E = 1
>>> sigma0 = L/50
>>> k0 = 50
>>> psi[:, 0] = np.exp(-(x - L/8)**2/(2*sigma0**2)) * np.exp(1j*k0*x)
>>> psi[:, 0] /= np.linalg.norm(psi[:, 0])
...
>>> (hv.Curve(np.absolute(psi[:, 0]), label='Wave function') *
>>>  hv.Curve(V, label='Potential'))
...
>>> # Trotter Suzuki
... for i in range(nt-1):
...     psi_p = np.fft.fft(np.exp(-1j*h*V/h_bar)*psi[:, i])
...     psi[:, i+1] = np.fft.ifft(np.exp(-1j*h*k**2/2/mass)*psi_p)
...
>>> %output holomap='scrubber'
>>> %output max_frames=100000
>>> hv.HoloMap([(
...             i*h, hv.Curve(np.real(psi[:, i])**2)) # * hv.Curve(V))
...             for i in range(nt)], kdims = ["Time"])
>>> #hv.HoloMap(np.absolute(psi), kdims=['Time'], label='Wave function')
b':HoloMap   [Time]\n   :Curve   [x]   (y)'
```

```python

```

```python

```

```python

```
