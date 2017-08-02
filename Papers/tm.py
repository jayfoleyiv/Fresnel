import cmath as cmath
import numpy as np
import numpy.linalg as npla

c = 299792458
theta = 45 * (cmath.pi) / 180
lpol = "s"
omega = 500e-9
numLayers = 3


def pmatrix(d, n, theta, omega,):
    kzl = cmath.sqrt(((n * omega / c) ** 2) - ((n * cmath.sin(theta) * omega / c) ** 2))
    print(kzl)
    phi = kzl * d
    P = np.array([[cmath.exp(-1j * phi), 0], [0, cmath.exp(1j * phi)]], dtype=complex)
    return P


def inversematrix(a, b, c, d):
    k = 1 / (a * d - b * c)
    matrix = k * np.array([d, -b], [-c, a])


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


def dmatrix(n, angle, lpol, layer):
    newTheta = t(angle, 1, layer)
    d = np.array  # does this do anything?
    if (lpol == "s"):
        d = np.array([[1 + 0j, 1 + 0j], [n * cmath.cos(newTheta), -n * cmath.cos(newTheta)]], dtype=complex)
    elif (lpol == "p"):
        d = np.array([[cmath.cos(newTheta), cmath.cos(newTheta)], [n, -n]], dtype=complex)
    return d


def getD(layer):
    if (layer == 1):
        return 0
    elif (layer == 2):
        return 0.5
    elif (layer == 3):
        return 0
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


def product(numLayers, lpol, theta, omega):
    sumproduct = np.identity(2, complex)
    for x in range(2, numLayers):
        d = getD(x)
        n = getN(x)
        dM = dmatrix(n, theta, lpol, x)
        pM = pmatrix(d, getN(1), theta, omega)
        print(pM)
        inverseDM = npla.inv(dM)  # check if the OOP works
        newMatrix = np.matmul(dM, np.matmul(pM, inverseDM))
        sumproduct *= newMatrix

    return sumproduct


def main():
    dM = dmatrix(getN(1), theta, lpol, 1)
    inverseDM = npla.inv(dM)
    sumproduct = product(numLayers, lpol, theta, omega)
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
