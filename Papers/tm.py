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


def matProduct(a , b):
    #c = np.matrix([[0 + 0j, 0 + 0j], [0 + 0j, 0 + 0j]])
    c = np.zeros((2,2), complex)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                c[i][j] += a[i][k] * b[k][j]
def product(numLayers, lpol, theta, omega):
    sumproduct = np.identity(2, complex)
    for x in range(2, numLayers):
        d = getD(x)
        n = getN(x)
        dM = dmatrix(n, theta, lpol, x)
        pM = pmatrix(d, n, theta, omega)
        inverseDM = npla.inv(dM)
        print()
        firstMatrix = np.dot(pM, inverseDM)
        print(firstMatrix)
        print(pM)
        newMatrix = matProduct(dM, firstMatrix)
        print(3)
        print(newMatrix)
        print(3)
        sumproduct *= newMatrix

    return sumproduct


def main():
    dM = dmatrix(getN(1), theta, lpol, 1)
    inverseDM = npla.inv(dM)
    sumproduct = product(numLayers, lpol, theta, omega)
    print("k")
    print(sumproduct)
    dL = getD(numLayers)
    nL = getN(numLayers)
    LDM = dmatrix(nL, theta, lpol, 3)
    Mmatrix = np.matmul(inverseDM, np.matmul(sumproduct, LDM))
    r = Mmatrix[1, 0] / Mmatrix[0, 0]
    R = r * np.conjugate(r)
    t = 1/Mmatrix[0,0]
    #mannually input N
    #finalAngle = t(theta, 1, numLayers)
    T = (t**2) * cmath.cos(theta)/(cmath.cos(theta))
    A = 1 - T - R
    print("s")
    print(R)
    print(T)
    print(A)
    return R

main()
