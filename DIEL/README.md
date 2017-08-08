Contains data files with dielectric function data as a function of wavelength for real materials.
Note refractive index (n) is related to dielectric function (eps) by the following:

eps = n^2
n = sqrt(eps)

both quantities are complex for materials that absorb light.
i.e. 
Re(eps) = Re(n)*Re(n) - Im(n)*Im(n)
Im(eps) = 2*Re(n)*Im(n)
 
For metals, which are strongly absorbing at all visible wavelengths, Im(n) > Re(n) and so 
Re(eps) tends to be negative for metals.  
