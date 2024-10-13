import math
import matplotlib.pyplot as plt
import numpy as np

media =0
mediaQuadrados =0
v =[]

f = open("dadosPendulo.txt", "r")
for i in f.readlines():
    i = float( i.strip())
    v.append(i)
    media = media +i
    mediaQuadrados = mediaQuadrados + i*i
 
f.close()    
    

media = media/100
mediaQuadrados = mediaQuadrados/100
sigmaQ = mediaQuadrados - media*media    

#Funcao Densidade de Probabilidade gaussiana
def gaussiana(x):
    return 1 / ((math.sqrt(2 * math.pi*sigmaQ))) * np.exp(-0.5 * ((x - media)**2)/sigmaQ)

x = np.linspace(media-3.5*math.sqrt(sigmaQ), media+3.5*math.sqrt(sigmaQ), 100)
plt.plot(x,gaussiana(x), color = "black", label = "gaussiana")


#plotagem histograma
nBins = math.ceil(1 + math.log(100,2))
plt.hist(v, bins= nBins, density = True, edgecolor = "k", color = "white",label="frequência exp.")


plt.xlabel("10T(s)")
plt.ylabel("Frequência")
plt.legend(loc = "upper left")

plt.show()

