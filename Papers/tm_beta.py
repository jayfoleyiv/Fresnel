import cmath as cmath
import numpy as np
import numpy.linalg as npla

c = 299792458
theta = 45 * (cmath.pi) / 180
lpol = "s"
omega = c * (cmath.pi * 2) / (500 * 10**-9)
numLayers = 3


def getD(layer):
    if (layer == 1):
        return 0 * 10**-6
    elif (layer == 2):
        return 0.5 * 10**-6
    elif (layer == 3):
        return 0 * 10**-6
    else:
        return "error"


def getN(layer):
    if (layer == 1):
        return 1 + 0j
    elif (layer == 2):
        return 0.8 + 2.5j
    elif (layer == 3):
        return 1 + 0j
    else:
        return "error"

def t(angle, ilayer, layer):
    # ilayer should always be 1
    if (ilayer == layer):
        return angle
    else:
        n1 = getN(ilayer)
        n2 = getN(ilayer + 1)
        sintheta = n1 * (cmath.sin(angle)/n2)
        angle = cmath.asin(sintheta)
        return t(angle, ilayer + 1, layer)

def pmatrix(d, n, theta, omega,):
    initN = getN(1)
    kz = cmath.sqrt(((n * omega / c) ** 2) - ((initN * cmath.sin(theta) * omega / c) ** 2))
    phi = kz * d
    P = np.array([[cmath.exp(-1j * phi), 0], [0, cmath.exp(1j * phi)]], dtype=complex)
    return P



def dmatrix(n, angle, lpol, layer):
    newTheta = t(angle, 1, layer)
    d = np.array  # does this do anything?
    if (lpol == "s"):
        d = np.array([[1 + 0j, 1 + 0j], [n * cmath.cos(newTheta), -n * cmath.cos(newTheta)]], dtype=complex)
    elif (lpol == "p"):
        d = np.array([[cmath.cos(newTheta), cmath.cos(newTheta)], [n, -n]], dtype=complex)
    return d


## loop based matrix multiplication
def matProduct(a , b):
    c = np.zeros((2,2), complex)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                c[i][j] += a[i][k] * b[k][j]

def product(numLayers, lpol, theta, omega):

    ## create complex arrays for outputs of matmult steps
    ## will complete product from right to left:
    ## lhs begins as identity
    ## step 1:  lhs * D = t1
    ## step 2:  t1  * P = t2
    ## step 3:  t2  * Di = t3
    ## step 4:  lhs = t3
    ## step 5:  repeat from step 1 until to nlayers

    lhs = np.identity(2,complex)
    t1  = np.zeros((2,2),complex)
    t2  = np.zeros((2,2),complex)
    t3  = np.zeros((2,2),complex)
    sumproduct = np.zeros((2,2), complex)
    for x in range(2, numLayers):
        ## build the base arrays for the current layer
        d = getD(x)
        n = getN(x)
        dM = dmatrix(n, theta, lpol, x)
        pM = pmatrix(d, n, theta, omega)
        inverseDM = npla.inv(dM)

        ## continue multiplication of previous product by 
        ## arrays of current layer
        np.dot(lhs, dM, t1)
        np.dot(t1, pM, t2)
        np.dot(t2, inverseDM, t3)
       
        ## save running product as lhs for next iteration
        lhs = t3.copy()
        #np.dot(pM, inverseDM, firstMatrix)
        #np.dot(dM,firstMatrix,sumproduct)

    ## copy final product to sumproduct    
    sumproduct = lhs.copy()
    return sumproduct


def main():
    ## create complex arrays for outputs of matmult steps
    temp = np.zeros((2,2),complex)
    Mmatrix = np.zeros((2,2),complex)
    # create initial D matrix and invert it
    dM = dmatrix(getN(1), theta, lpol, 1)
    inverseDM = npla.inv(dM)
 
    # create product of D P D^-1 for all intermediate layers
    sumproduct = product(numLayers, lpol, theta, omega)
    
    # create D for final layer
    dL = getD(numLayers)
    nL = getN(numLayers)
    LDM = dmatrix(nL, theta, lpol, 3)

    # This is SP = PRODUCT_l (D_l P D_l^-1) * D_Nlayers
    np.dot(sumproduct, LDM, temp)
     
    # This is M = D_1 * SP
    np.dot(inverseDM, temp, Mmatrix)

    # Now get the r, t amplitudes along with R and T
    r = Mmatrix[1, 0] / Mmatrix[0, 0]
    R = r * np.conjugate(r)
    t = 1/Mmatrix[0,0]

    #mannually input N
    #finalAngle = t(theta, 1, numLayers)
    T = (t**2) * cmath.cos(theta)/(cmath.cos(theta))
    A = 1 - T - R
    print("M Matix")
    print(Mmatrix)
    print(" Reflection ")
    print(R)
    print(" Transmission ")
    print(T)
    print(" Absorption ")
    print(A)
    return R

main()

