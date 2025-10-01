import csv
import matplotlib.pyplot as plt
import math



lOna = 1.5406           #Angstrom

angle = []          #angle = 2 theta (degree)
intensity = []

with open('entrega1/files/archivo.csv') as file:
    reader = csv.reader(file)
    next(reader)
    for row in reader:
        x, y = row
        angle.append(float(x))
        intensity.append(float(y))

plt.style.use('classic')
plt.figure(figsize=(6,6), facecolor='white')
plt.plot(angle, intensity,'o', label='Experimental', color='red')
#plt.xlim( , )
#plt.ylim( , )
plt.xlabel(r'$2\theta$ [$\degree$]')
plt.ylabel(r'$I_{rel}$')
plt.grid(True)
plt.legend(loc = 'lower right')
plt.savefig('entrega1/informe/images/grafic_experimental.png', dpi=300, bbox_inches='tight')
plt.show()
