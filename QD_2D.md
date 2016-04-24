```python
>>> import numpy as np
>>> import scipy as sp
>>> import holoviews as hv
>>> hv.notebook_extension()
<IPython.core.display.HTML object>
```

```python
>>> L = 2
>>> nx = 256
>>> nt = 200
>>> h_bar = 1
>>> mass= 1/2
>>> a = L/(nx - 1) # space step
>>> h = 1e-2 # time step
>>> #x = np.linspace(0, L, nx)
... dk = 2*np.pi/L
>>> k0 = 0
>>> k = np.fft.fftfreq(nx, a)
>>> kx, ky = np.meshgrid(k, k)
>>> kr = kx**2+ky**2
>>> print(kr.shape)
>>> X, Y = np.mgrid[0:L+a:a, 0:L+a:a]
>>> print(X)
...
>>> # Potential Well
... V0 = 1500
>>> Vwidth = 0.01
>>> Vcenter = 0.6*L
>>> V = np.piecewise(x, [x < Vcenter - Vwidth/2, x >= Vcenter - Vwidth/2, x > Vcenter + Vwidth/2], [0, V0, 0])
>>> V = np.zeros((nx, nx))
...
>>> # Initial wave packet
... psi = np.zeros((nx, nx, nt), dtype='cfloat')
>>> sigma0 = 40
>>> psi[:, :, 0] = np.exp(-(X - L/2)**2/(2*sigma0**2)) * np.exp(-(Y - L/2)**2/(2*sigma0**2)) * np.exp(1j*k0*X)
>>> psi[:, :, 0] /= np.linalg.norm(psi[:, :, 0])
>>> (hv.Image(np.absolute(psi[:, :, 0])**2, label='Wave function')) #*
>>>  #hv.Image(V, label='Potential'))
(256, 256)
[[ 0.          0.          0.         ...,  0.          0.          0.        ]
 [ 0.00784314  0.00784314  0.00784314 ...,  0.00784314  0.00784314
   0.00784314]
 [ 0.01568627  0.01568627  0.01568627 ...,  0.01568627  0.01568627
   0.01568627]
 ..., 
 [ 1.98431373  1.98431373  1.98431373 ...,  1.98431373  1.98431373
   1.98431373]
 [ 1.99215686  1.99215686  1.99215686 ...,  1.99215686  1.99215686
   1.99215686]
 [ 2.          2.          2.         ...,  2.          2.          2.        ]]
b':Image   [x,y]   (z)'
```

## Trotter-Suzuki

```python
>>> # Trotter Suzuki
... for i in range(nt-1):
...     psi_p = np.fft.fft2(np.exp(-1j*h*V/h_bar)*psi[:, :, i])
...     psi[:, :, i+1] = np.fft.ifft2(np.exp(-1j*h*kr/2/mass)*psi_p)
```

```python
>>> %output holomap='scrubber'
>>> %output max_frames=100000
>>> hv.HoloMap([(
...             i*h, hv.Image(np.absolute(psi[:, :, i])**2)) # * hv.Curve(V))
...             for i in range(nt)], kdims = ["Time"])
>>> #hv.HoloMap(np.absolute(psi), kdims=['Time'], label='Wave function')
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
