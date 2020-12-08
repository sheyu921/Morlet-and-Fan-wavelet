"""
calculate the wavenumber from spatial data
input
matrix: data matrix (topography or gravity)
dx, dy: the interval in x and y direction

output
     K: 2D wavenumbers
kx, ky: wavenumbers in x and y direction
"""
def xy2kxky(matrix,dx,dy):
    import numpy as np
    [M,N] = matrix.shape
    kx = 2*np.pi*np.linspace(-np.floor(N/2),np.floor((N-1)/2),N)/N/dx
    ky = 2*np.pi*np.linspace(-np.floor(M/2),np.floor((M-1)/2),M)/M/dy
    [kx,ky] = np.meshgrid(kx,ky)
    K = np.sqrt(kx**2+ky**2)   
    return K, kx, ky
    
"""
calculate the Morlet wavelet
input
kx, ky: wavenumbers in x and y direction
     s: scale of the wavelet
 theta: the degree from the east axis direction
    k0: the central wavenumber

output
psik: spectral morlet wavelet
"""
def spectral_morlet(kx,ky,s,theta,k0):
    import numpy as np
    theta = np.radians(theta)
    a = (kx*s-k0*np.cos(theta))**2
    b = (ky*s-k0*np.sin(theta))**2
    if k0 > 5:
        psik = np.exp(-(a+b)/2)*s
    else:
        psik = (np.exp(-(a+b)/2)-np.exp(-((kx*s)**2+(ky*s)**2+k0**2)/2))*s
    return psik

"""
calculate the fan wavelet
input
  kx, ky: wavenumbers in x and y direction
       s: scale of the wavelet
      k0: the central wavenumber

output
psik_fan: spectral fan wavelet
"""
def fan(kx, ky, s, k0, sp):
    import numpy as np
    dtheta = np.rad2deg(2*np.sqrt(-2*np.log(0.75))/k0)
    Ntheta = int(np.round(180/dtheta))
    thetac = 90
    theta0 = thetac-(Ntheta-1)*dtheta/2
    psik_fan = np.zeros_like(kx)
    for i in range(Ntheta):
        theta = i*dtheta+theta0;
        psik_tmp = spectral_morlet(kx, ky, s, theta, k0);
        psik_fan = psik_fan + psik_tmp
    return psik_fan
