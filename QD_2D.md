## 2D Quantum Dynamics

```python
>>> import numpy as np
>>> import scipy as sp
>>> import holoviews as hv
>>> hv.notebook_extension()
<IPython.core.display.HTML object>
```

```python
>>> L = 2 # x domain from 0 to L
>>> nx = 251
>>> nt = 2000
>>> h_bar = 1
>>> mass= 1/2
>>> a = L/(nx - 1) # space step
>>> h = 1e-2 # time step
>>> x = np.linspace(0, L, nx)
>>> dk = 2*np.pi/L
>>> k0 = 50
>>> k = np.linspace(-np.pi/L, np.pi/L, num=nx)
...
>>> kx, ky = np.meshgrid(k, k)
>>> #kx, ky = np.mgrid[-np.pi/L:np.pi/L:dk, -np.pi/L:np.pi/L:dk]
... k2 = kx**2 + ky**2
>>> X, Y = np.mgrid[0:L+a:a, 0:L+a:a]
...
>>> # Potential Well
... V0 = 4200
>>> V = np.zeros((nx, nx))
>>> V[100:105, 0:96] = V0
>>> V[100:105, 120:132] = V0
>>> V[100:105, 156:nx] = V0
...
>>> # Initial wave packet
... psi = np.zeros((nx, nx, nt), dtype='cfloat')
>>> sigma0 = 0.1
>>> psi[:, :, 0] = np.exp(-(X - L/2)**2/(2*sigma0**2)) * np.exp(-(Y - L/2)**2/(2*sigma0**2)) * np.exp(1j*k0*X)
>>> psi[:, :, 0] /= np.linalg.norm(psi[:, :, 0])
>>> %opts Image norm{+axiswise}
>>> (hv.Image(np.absolute(psi[:, :, 0])**2, label='Wave function') + hv.Image(V, label='Potential'))
:Layout
   .Image.Wave_function :Image   [x,y]   (z)
   .Image.Potential     :Image   [x,y]   (z)
```

## Trotter-Suzuki

```python
>>> # Trotter Suzuki
... for i in range(nt-1):
...     psi_p = np.fft.fft2(np.exp(-1j*h*V/h_bar)*psi[:, :, i])
...     psi[:, :, i+1] = np.fft.ifft2(np.exp(-1j*h*k2/2/mass)*psi_p)
```

```python
>>> %output holomap='scrubber'
>>> %output max_frames=100000
>>> hv.HoloMap([(
...             i*h, hv.Image(np.absolute(psi[:, :, i])**2))
...             for i in range(nt)], kdims = ["Time"])
b':HoloMap   [Time]\n   :Image   [x,y]   (z)'
```

```python

```

```python

```

```python

```

```python

```
