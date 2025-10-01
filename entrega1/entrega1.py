import csv
import matplotlib.pyplot as plt
import math



lOna = 1.5406           #Angstrom

angle = []          #angle = 2 theta (degree)
intensity = []

with open("archivo.csv") as file:
    reader = csv.reader(file)
    next(reader)
    for row in reader:
        x, y = row
        angle.append(float(x))
        intensity.append(float(y))