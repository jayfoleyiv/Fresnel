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

#### 
Added on August 24th:  Material constants for other high-temperature metals and alloys.
All data sets are structured as follows:  Col1 -> wavelength in meters, Col2 -> Re(eps), Col3 -> Im(eps)
- Cr_Palik.txt -> Chromium (melting point 2180 K)
- Pd_Palik.txt -> Palladium (melting point 1828 K)
- Pt_Palik.txt -> Platinum (melting point 2041 K)
- Rh_Palik.txt -> Rhodium (melting point 2236 K)
- Ta_CRC.txt   -> Tantalum (melting point 3293 K)
- TiN_Palik.txt -> Titanium nitride (melting point 3203 K)
- Ti_Palik.txt  -> Titanium (melting point 1941 K) 
- V_CRC.txt -> Vanadium (melting point 2183 K)

