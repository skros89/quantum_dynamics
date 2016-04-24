## 2D Quantum Dynamics

```python
>>> import numpy as np
>>> import scipy as sp
>>> import holoviews as hv
>>> hv.notebook_extension()
<IPython.core.display.HTML object>
```

```python
>>> L = 2 # x and y domain from 0 to L
>>> nx = 251
>>> nt = 3000
>>> a = L/(nx - 1) # space step
>>> h = 1e-2 # time step
>>> x = np.linspace(0, L, nx)
>>> #dk = 2*np.pi/L
... k0 = 50
>>> k = np.linspace(-np.pi/L, np.pi/L, num=nx)
...
>>> kx, ky = np.meshgrid(k, k)
>>> #kx, ky = np.mgrid[-np.pi/L:np.pi/L:dk, -np.pi/L:np.pi/L:dk]
... k2 = kx**2 + ky**2
>>> X, Y = np.mgrid[0:L+a:a, 0:L+a:a]
...
>>> # Double slit potential
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
>>> (hv.Image(np.absolute(psi[:, :, 0])**2, label='Probability density of the Wave function') + hv.Image(V, label='Potential'))
:Layout
   .Image.Wave_function :Image   [x,y]   (z)
   .Image.Potential     :Image   [x,y]   (z)
```

```python
>>> # Trotter Suzuki
... for i in range(nt-1):
...     psi_p = np.fft.fft2(np.exp(-1j*h*V)*psi[:, :, i])
...     psi[:, :, i+1] = np.fft.ifft2(np.exp(-1j*h*k2)*psi_p)
```

```python
>>> %%output filename="2d_double_slit" fig="png" holomap="gif"
... %output holomap='scrubber'
... %output max_frames=100000
... hv.HoloMap([(i*h, hv.Image(np.absolute(psi[:, :, i])**2))for i in range(nt)], kdims = ["Time"], label='Probability density of the Wave function'))
```
