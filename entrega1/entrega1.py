import csv
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np


################################################################################


class RegLin:
    def __init__(self, _x, _y):

        self.NombreMesures = len(_x)
        sumX = 0
        sumY = 0
        sumX2= 0
        sumY2= 0
        sumXY= 0
        for n in range(self.NombreMesures):
            sumX = sumX + _x[n]
            sumY = sumY + _y[n]
            sumX2 += _x[n] * _x[n]
            sumY2 += _y[n] * _y[n]
            sumXY += _x[n] * _y[n]

        deltaX = (self.NombreMesures*sumX2)-(sumX*sumX)
        deltaY = (self.NombreMesures*sumY2)-(sumY*sumY)

        self.pendent = ((self.NombreMesures*sumXY)-(sumX*sumY))/deltaX
        self.ordenadaOrigen = ((sumY*sumX2)-(sumX*sumXY))/deltaX

        aux = 0
        for n in range(self.NombreMesures):
            aux += pow(_y[n] - self.ordenadaOrigen - (self.pendent*_x[n]), 2)
        o2 = (1 /(self.NombreMesures-2))*aux

        self.errorPendent = np.sqrt((o2/deltaX) * self.NombreMesures)
        self.errorOrdenadaOrigen = np.sqrt((o2/deltaX) * sumX2)

        self.R2 = (pow((self.NombreMesures*sumXY) - (sumX*sumY),2))/(deltaX * deltaY)

    def Calcular(self, x):
        return (self.pendent * x) + self.ordenadaOrigen

    def ShowInfo(self):
        print("Pendent: ", self.pendent, " ± ", self.errorPendent, "\nOrdenada a origen: ", self.ordenadaOrigen, " ± ", self.errorOrdenadaOrigen, "\nR^2: ", self.R2)


################################################################################


waveLength = 1.5406           #Angstrom

angle = []          #angle = 2 theta (degree)
intensity = []

with open('entrega1/files/19_diffrac.csv') as file:
    reader = csv.reader(file)
    next(reader)
    for row in reader:
        x, y = row
        angle.append(float(x))
        intensity.append(float(y))

plt.style.use('classic')
plt.figure(figsize=(6,6), facecolor='white')
plt.plot(angle, intensity, label='Experimental', color='red')
#plt.xlim( , )
plt.ylim(-50,350)
plt.xlabel(r'$2\theta$ [$\degree$]')
plt.ylabel(r'$I_{relative}$ [counts]')
plt.grid(True)
plt.legend(loc = 'lower right')
plt.savefig('entrega1/informe/images/grafic_experimental.png', dpi=300, bbox_inches='tight')
#plt.show()

index, _ = find_peaks(intensity, height = 50)
anglePeaks = [angle[i] for i in index]
intensityPeaks = [intensity[i] for i in index]

plt.style.use('classic')
plt.figure(figsize=(6,6), facecolor='white')
plt.plot(angle, intensity, label='Experimental', color='red')
plt.plot(anglePeaks, intensityPeaks, 'o', label='Peaks', color='blue')
#plt.xlim( , )
plt.ylim(-50,350)
plt.xlabel(r'$2\theta$ [$\degree$]')
plt.ylabel(r'$I_{relative}$ [counts]')
plt.grid(True)
#plt.legend(loc = 'lower right')
plt.savefig('entrega1/informe/images/grafic_experimental_peaks.png', dpi=300, bbox_inches='tight')
#plt.show()


################################################################################


dExperimental = [waveLength / (2 * np.sin( np.radians( twoTheta / 2 ) )) for twoTheta in anglePeaks]
factorD = [1 / (d ** 2) for d in dExperimental]


allVectors = [[h, k, l, mod] for l in range(10) for k in range(10) for h in range(10) for mod in [h**2 + k**2 + l**2]
           if (h,k,l) != (0,0,0)]
sortedVectors = sorted(allVectors, key=lambda x: x[3])
vectorsNoRepeat = []
for vec in sortedVectors:
    if vec[3] not in [v[3] for v in vectorsNoRepeat]:
        vectorsNoRepeat.append(vec)
vectors = vectorsNoRepeat[:30]


def StructureFactorSC(G):
    h, k, l, mod = G
    c_0, c_1, c_2, c_3 = [0,0,0], [1,0,0], [0,1,0], [0,0,1]
    c = [c_0, c_1, c_2, c_3]
    S = 0
    for c_i in c:
        S += np.exp(2j * np.pi * (h*c_i[0] + k*c_i[1] + l*c_i[2]))
    return S
def StructureFactorBCC(G):
    h, k, l, mod = G
    c_0, c_1 = [0,0,0], [0.5, 0.5, 0.5]
    c = [c_0, c_1]
    S = 0
    for c_i in c:
        S += np.exp(2j * np.pi * (h*c_i[0] + k*c_i[1] + l*c_i[2]))
    return S
def StructureFactorFCC(G):
    h, k, l, mod = G
    c_0, c_1, c_2, c_3 = [0,0,0], [0.5,0.5,0], [0.5,0,0.5], [0,0.5,0.5]
    c = [c_0, c_1, c_2, c_3]
    S = 0
    for c_i in c:
        S += np.exp(2j * np.pi * (h*c_i[0] + k*c_i[1] + l*c_i[2]))
    return S
def StructureFactorDiamond(G):
    h, k, l, mod = G
    c_00, c_10, c_20, c_30 = [0,0,0], [0.5,0.5,0], [0.5,0,0.5], [0,0.5,0.5]
    c_01, c_11, c_21, c_31 = [0.25,0.25,0.25], [0.75,0.75,0.25], [0.75,0.25,0.75], [0.25,0.75,0.75]
    c = [c_00, c_10, c_20, c_30, c_01, c_11, c_21, c_31]
    S = 0
    for c_i in c:
        S += np.exp(2j * np.pi * (h*c_i[0] + k*c_i[1] + l*c_i[2]))
    return S


modSC = [v[3] for v in vectors if StructureFactorSC(v) != 0]
modBCC = [v[3] for v in vectors if StructureFactorBCC(v) != 0]
modFCC = [v[3] for v in vectors if StructureFactorFCC(v) != 0]
modDiamond = [v[3] for v in vectors if StructureFactorDiamond(v) != 0]
mods = [modSC[:len(factorD)], modBCC[:len(factorD)], modFCC[:len(factorD)], modDiamond[:len(factorD)]]


################################################################################


for mod, name in zip(mods, ['SC', 'BCC', 'FCC', 'Diamond']):

    regLin = RegLin(factorD, mod)
    x_min, x_max = 0.1, 0.7

    plt.style.use('classic')
    plt.figure(figsize=(6,6), facecolor='white')
    plt.plot(factorD, mod, 'o', label=r'$\frac{1}{d^2}(h^2 + k^2 + l^2)$', color='red')
    plt.plot([x_min, x_max], [regLin.Calcular(x_min), regLin.Calcular(x_max)], label=f'$R^2 = {regLin.R2:.3f}$', color='blue')
    #plt.xlim( , )
    plt.ylim(0.5,6.5)
    plt.xlabel(r'$2\theta$ [$\degree$]')
    plt.ylabel(r'$I_{relative}$ [counts]')
    plt.grid(True)
    plt.legend(loc = 'lower right')
    plt.savefig(f'entrega1/informe/images/{name}.png', dpi=300, bbox_inches='tight')
    #plt.show()