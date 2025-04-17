import pandas
import numpy
import scipy

filename = input("File name: ")
guessFile = pandas.read_csv(filename)
temperatures = guessFile["T (K)"].tolist()
pressures = guessFile["P (Mpa)"].tolist()

liqWaterTemps = []
iceWaterTemps = []

liqWaterPressures = []
iceWaterPressures = []

for i in range(len(temperatures)):
    if temperatures[i] >= 273.15:
        liqWaterTemps.append(temperatures[i])
        liqWaterPressures.append(pressures[i]*1E6)
    else:
        iceWaterTemps.append(temperatures[i])
        iceWaterPressures.append(pressures[i]*1E6)
        
def generateParameters(T, P):
    def model(T, A, B):
        return A*numpy.exp(B*T)
    
    initialGuess = [1E-10, 0.1]
    
    A, B = scipy.optimize.curve_fit(model, numpy.array(T).flatten(), numpy.array(P).flatten(), p0=initialGuess, maxfev=5000)[0]
    return A, B

try:
    A1, B1 = generateParameters(iceWaterTemps, iceWaterPressures)
    print("A1: " + str(A1))
    print("B1: " + str(B1))
except:
    print("Ice phase parameters could not be generated")
A2, B2 = generateParameters(liqWaterTemps, liqWaterPressures)
print("A2: " + str(A2))
print("B2: " + str(B2))

