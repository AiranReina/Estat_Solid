import csv
import matplotlib.pyplot as plt
import numpy as np

print('')
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

def mfp(K, c_V):
    return [3 * K[i] / (c * c_V[i]) for i in range(len(K))] #[m]

def conductivity(L,Temp,c_v,T):
    for i,t in enumerate(T):
        if t == Temp:
            j = i
    return (1/3) * c * L * c_v[j] #[W/mK]

################################################################################

c = 6400 #[m/s]
kB = 1.38e-23 #[J/K]
h = 6.63e-34 #[J·s]
Na = 6.022e23 #[atoms/mol]

mSi = 28 #[uma]
rhoSi = 2330 #[kg/m3]
nSi = (rhoSi/mSi) * 1000 * Na  #[atoms/m3]

T = [3,4,5,6,7,8,9,10,20,30,50,70,90,100,300,500,700,900,1100] #[K]
K = [129,271,494,816,1220,1675,2150,2623,5446,4968,2797,1618,1060,845,136,70,48,36,30.5] #[W/mK]
c_v = [0.0011,0.0027,0.0053,0.0092,0.0145,0.0217,0.0309,0.0424,0.3389,1.1354,4.5175,
       8.8614,12.6303,14.1650,23.2528,24.3145,24.6191,24.7461,24.8107] #[J·K/mol]
c_v_vol = [cv * nSi / Na for cv in c_v] #[J·K/m3]

################################################################################

plt.style.use('classic')
plt.figure(figsize=(6,6), facecolor='white')
plt.plot(T, K , color='blue')
plt.plot(T, K , 'o', color='blue', label=r'Conductivitat tèrmica $\lambda$')
plt.xlabel(r'$T$ [K]')
plt.ylabel(r'$\lambda$ [W/mK]')
plt.grid(True)
plt.legend(loc = 'upper right')
#plt.title()
plt.savefig(f'entrega2/informe/images/lambda_T.png', dpi=300, bbox_inches='tight')
#plt.show()

plt.style.use('classic')
plt.figure(figsize=(6,6), facecolor='white')
plt.plot(T, c_v, color='red')
plt.plot(T, c_v, 'o', color='red', label=r'Capacitat calorífica $c_v$')
plt.xlabel(r'$T$ [K]')
plt.ylabel(r'$c_v$ [J·K/mol]')
plt.grid(True)
plt.legend(loc = 'lower right')
#plt.title()
plt.savefig(f'entrega2/informe/images/cv_T.png', dpi=300, bbox_inches='tight')
#plt.show()

################################################################################

lim = 8
T3 = [t**3 for t in T[0:lim]]

regLin = RegLin(T3, c_v)
x_min, x_max = 0, T3[-1]
theta_D = np.cbrt(((12 * np.pi**4)/5) * kB * Na / regLin.pendent)

plt.style.use('classic')
plt.figure(figsize=(6,6), facecolor='white')
plt.plot([x_min, x_max], [regLin.Calcular(x_min), regLin.Calcular(x_max)], label=f'$R^2 = {regLin.R2:.3f}$', color='black')
plt.plot(T3, c_v[0:lim], 'o', color='red', label=r'Capacitat calorífica $c_v$')
plt.xlabel(r'$T^3$ [K$^3$]')
plt.ylabel(r'$c_v$ [J·K/mol]')
plt.grid(True)
plt.legend(loc = 'lower right')
#plt.title()
plt.savefig(f'entrega2/informe/images/cv_T3.png', dpi=300, bbox_inches='tight')
#plt.show()

print(f'Theta_D = {theta_D} K')
#print(f'Pendent = {regLin.pendent} ± {regLin.errorPendent} J·K/mol·K^3')

################################################################################

l = mfp(K, c_v_vol) #[m]
log_l = [np.log(l_i) for l_i in l]
mfp_mm = [l_i * 1e3 for l_i in l] #[mm]
#T_inv = [ti**(-1.595) for ti in T]
log_T = [np.log(Ti) for Ti in T]

plt.style.use('classic')
plt.figure(figsize=(6,6), facecolor='white')
plt.plot(T, mfp_mm, color='green')
plt.plot(T, mfp_mm, 'o', color='green', label=r'Camí lliure mitjà $l$')
plt.xlabel(r'$T$ [K]')
plt.ylabel(r'$l$ [mm]')
plt.grid(True)
plt.legend(loc = 'upper right')
#plt.title()
plt.savefig(f'entrega2/informe/images/mfp_T.png', dpi=300, bbox_inches='tight')
#plt.show()

lim2 = 14
plt.style.use('classic')
plt.figure(figsize=(6,6), facecolor='white')
plt.plot(T[2:lim2], mfp_mm[2:lim2], color='green')
plt.plot(T[0:lim2], mfp_mm[0:lim2], 'o', color='green', label=r'Camí lliure mitjà $l$')
plt.plot(T[:2], mfp_mm[:2], 'o', color='red')
plt.xlabel(r'$T$ [K]')
plt.ylabel(r'$l$ [mm]')
plt.grid(True)
plt.legend(loc = 'upper right')
#plt.title()
plt.savefig(f'entrega2/informe/images/mfp_T_zoomed.png', dpi=300, bbox_inches='tight')
#plt.show()

lim3 = 15
regLin = RegLin(log_T[lim3:], log_l[lim3:])
x_min, x_max = log_T[lim3], log_T[-1]
plt.style.use('classic')
plt.figure(figsize=(6,6), facecolor='white')
plt.plot([x_min, x_max], [regLin.Calcular(x_min), regLin.Calcular(x_max)], label=f'$R^2 = {regLin.R2:.3f}$', color='black')
plt.plot(log_T[lim3:], log_l[lim3:], 'o', color='green', label=r'Camí lliure mitjà $l$')
plt.xlabel(r'$log(T)$ [K]')
plt.ylabel(r'$log(l)$ [mm]')
plt.grid(True)
plt.legend(loc = 'upper right')
plt.text(6.3, -18.9, f'Pendent = {regLin.pendent:.3f} ± {regLin.errorPendent:.3f}', fontsize=10)
#plt.title()
plt.savefig(f'entrega2/informe/images/mfp_T_logarithmic.png', dpi=300, bbox_inches='tight')
#plt.show()
print(f'Pendent = {regLin.pendent} ± {regLin.errorPendent}')

################################################################################

Conductivitats = []
Temperatures = [10,100,500]
Mides = [1e-6, 100e-9, 10e-9]

for L in Mides:
    Conductivitats.append([])
    for Temp in Temperatures:
        Conductivitats[-1].append(conductivity(L,Temp,c_v_vol,T))
print('')
print(Conductivitats)
print('')

################################################################################

K_total = []
for i,Ti in enumerate(T):
    if Ti == 10 or Ti == 100 or Ti == 500:
        print(mfp_mm[i] * 1000, r'\mu m at', Ti, 'K')
        K_total.append([])
        for L in Mides:
            K_boundary = (1/3) * c * L * c_v_vol[i]
            K_bulk = K[i]
            K_total[-1].append(1 / ((1 / K_bulk) + (1 / K_boundary)))

print('')
print(K_total)

for L in Mides:
    mfp_L = []
    for li in l:
        mfp_L.append(1e9 / ((1 / li) + (1 / L)))
    plt.style.use('classic')
    plt.figure(figsize=(6,6), facecolor='white')
    plt.plot(T, mfp_L, color='green')
    plt.plot(T, mfp_L, 'o', color='green', label=r'Camí lliure mitjà $l$')
    plt.xlabel(r'$T$ [K]')
    plt.ylabel(r'$l$ [nm]')
    plt.grid(True)
    plt.legend(loc = 'upper right')
    #plt.title()
    plt.savefig(f'entrega2/informe/images/mfp_T_{L}.png', dpi=300, bbox_inches='tight')
    #plt.show()

    plt.style.use('classic')
    plt.figure(figsize=(6,6), facecolor='white')
    plt.plot(T, mfp_L, color='green')
    plt.plot(T, mfp_L, 'o', color='green', label=r'Camí lliure mitjà $l$')
    plt.xlabel(r'$T$ [K]')
    plt.ylabel(r'$l$ [nm]')
    plt.xlim(0,200)
    plt.grid(True)
    plt.legend(loc = 'lower left')
    #plt.title()
    plt.savefig(f'entrega2/informe/images/mfp_T_{L}_zoomed.png', dpi=300, bbox_inches='tight')
    #plt.show()