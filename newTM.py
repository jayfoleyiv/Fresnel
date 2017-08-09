import cmath as cmath
import numpy as np
import numpy.linalg as npla

c = 299792458
theta = 45 * (cmath.pi) / 180
lpol = "p"
#omega = c * (cmath.pi * 2) / (500 * 10 ** -9)
numLayers = 3


arr = np.loadtxt("/Users/jtsatsaros2018/Fresnel/FDTD_TEST/Fresnel/DIEL/W_Bulk_epsilon_vs_Wavelength_in_Meters.txt")
#ext = ("/Users/jtsatsaros2018/Fresnel/FDTD_TEST/Fresnel/DIEL/W_Bulk_epsilon_vs_Wavelength_in_Meters.txt")
#for line in (text):
    #linecount = 0
    #m = text.readline(linecount)
    #l = []
    #st = ""
   # print(m)
  #  for n in m:
        #if (m[n] != " "):
       #     st.join(m[n])
      #  else:
     #       i = int(st)
    #        l.append(i)
   #         st = ""
  #  arr.append(l)
 #   linecount += 1
#text.close()


def getD(layer):
    if (layer == 1):
        return 0 * 10 ** -6
    elif (layer == 2):
        return 10 * 10 ** -9
    elif (layer == 3):
        return 0 * 10 ** -6
    else:
        return "error"

def getN(diearray , layer, index):
    if (layer == 1):
        return 1 + 0j
    elif (layer == 2):
        return cmath.sqrt(int(diearray[index][1]) + int(diearray[index][2])*1j)
    elif (layer == 3):
        return 1 + 0j
    else:
        return "error"

def t(angle, ilayer, layer, a , z):
    # ilayer should always be 1
    if (ilayer == layer):
        return angle
    else:
        n1 = getN(a, ilayer, z)
        n2 = getN(a, ilayer + 1, z)
        sintheta = n1 * (cmath.sin(angle) / n2)
        angle = cmath.asin(sintheta)
        return t(angle, ilayer + 1, layer, a , z)


def pmatrix(d, n, angle, omega, a, z):
    initN = getN(a, 1, z)
    kz = cmath.sqrt(((n * omega / c) ** 2) - ((initN * cmath.sin(angle) * omega / c) ** 2))
    phi = kz * d
    P = np.array([[cmath.exp(-1j * phi), 0], [0, cmath.exp(1j * phi)]], dtype=complex)
    return P


def dmatrix(n, angle, lpol, layer, a, z):
    newTheta = t(angle, 1, layer, a, z)
    d = np.array  # does this do anything?
    if (lpol == "s"):
        d = np.array([[1 + 0j, 1 + 0j], [n * cmath.cos(newTheta), -n * cmath.cos(newTheta)]], dtype=complex)
    elif (lpol == "p"):
        d = np.array([[cmath.cos(newTheta), cmath.cos(newTheta)], [n, -n]], dtype=complex)
    return d


## loop based matrix multiplication
def matProduct(a, b):
    c = np.zeros((2, 2), complex)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                c[i][j] += a[i][k] * b[k][j]


def product(numLayers, lpol, angle, omega, a, z):
    ## create complex arrays for outputs of matmult steps
    ## will complete product from right to left:
    ## lhs begins as identity
    ## step 1:  lhs * D = t1
    ## step 2:  t1  * P = t2
    ## step 3:  t2  * Di = t3
    ## step 4:  lhs = t3
    ## step 5:  repeat from step 1 until to nlayers

    lhs = np.identity(2, complex)
    t1 = np.zeros((2, 2), complex)
    t2 = np.zeros((2, 2), complex)
    t3 = np.zeros((2, 2), complex)
    sumproduct = np.zeros((2, 2), complex)
    for x in range(2, numLayers):
        ## build the base arrays for the current layer
        d = getD(x)
        n = getN(a, x, z)
        dM = dmatrix(n, angle, lpol, x, a, z)
        pM = pmatrix(d, n, angle, omega, a, z)
        inverseDM = npla.inv(dM)

        ## continue multiplication of previous product by
        ## arrays of current layer
        np.dot(lhs, dM, t1)
        np.dot(t1, pM, t2)
        np.dot(t2, inverseDM, t3)

        ## save running product as lhs for next iteration
        lhs = t3.copy()
        # np.dot(pM, inverseDM, firstMatrix)
        # np.dot(dM,firstMatrix,sumproduct)

    ## copy final product to sumproduct
    sumproduct = lhs.copy()
    return sumproduct


def main():
    ar = np.array(["wavelength", "reflectance", "transmittance"])
    z = 0
    for y in range(len(arr)):
        wl = arr[y][0]
        om = c * (cmath.pi * 2) / (wl)

        ## create complex arrays for outputs of matmult steps
        temp = np.zeros((2, 2), complex)
        Mmatrix = np.zeros((2, 2), complex)
        # create initial D matrix and invert it
        dM = dmatrix(getN(arr, 1, z), theta, lpol, 1, arr, z)
        inverseDM = npla.inv(dM)

        # create product of D P D^-1 for all intermediate layers
        sumproduct = product(numLayers, lpol, theta, om, arr, z)
        # create D for final layer
        dL = getD(numLayers)
        nL = getN(arr, numLayers, z)
        LDM = dmatrix(nL, theta, lpol, 3, arr, z)

        # This is SP = PRODUCT_l (D_l P D_l^-1) * D_Nlayers
        np.dot(sumproduct, LDM, temp)

        # This is M = D_1 * SP
        np.dot(inverseDM, temp, Mmatrix)
        # Now get the r, t amplitudes along with R and T
        r = Mmatrix[1, 0] / Mmatrix[0, 0]
        R = r * np.conjugate(r)
        lt = 1 / Mmatrix[0, 0]

        finalAngle = t(theta, 1, numLayers, arr, z)



        T = lt*np.conj(lt) * (getN(arr, numLayers, z)) * cmath.cos(finalAngle) / ((cmath.cos(theta) * getN(arr, 1, z)))
        A = 1 - T - R
        m = np.array([wl, R, T])
        ar = np.vstack((ar,m))
        ar = np.real(ar)
        z +=1

        ans = [A, T, R]

  #  for k in range(len(ar)):
   #     for l in range(3):
    #        a = ar[k][l]
     #       a = np.real(a)
      #      print(a, end =" ")
       # print()

    text_file = open("/Users/jtsatsaros2018/Documents/test1", "a")
    for k in range(len(ar)):
        for l in range(3):
            a = ar[k][l]
            a = np.real_if_close(a)
            text_file.write(str(a))
        text_file.write("\n")
    text_file.close()

    return ans


print(main())

