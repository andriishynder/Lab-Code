# %%
import matplotlib.pyplot as plt
import numpy as np
from time import sleep

from pydaqmx_helper.dac import DAC
from pydaqmx_helper.adc import ADC
from pydaqmx_helper.digital_io import Digital_IO     

myDAC = DAC(0)

myADC = ADC()
myADC.addChannels([1])

myDigital_IO = Digital_IO('0','output')

# %%
def StepMotorOuter(CW):
    if CW:
        myDigital_IO.write(0b10011)
    else:
        myDigital_IO.write(0b10001)
        
    sleep(0.06)
    myDigital_IO.write(0b00000)

    return

def StepMotorInner(CW):
    if CW:
        myDigital_IO.write(0b01011)
    else:
        myDigital_IO.write(0b01001)
   
    sleep(0.1)
    myDigital_IO.write(0b00000) 

    return

# %%
data = [] 
data2 = [] 
data3 = [] 
myDigital_IO.write(0b00000)

DegToStepOuter = 25
DegToStepInner = 25
CW = False;
#CW = True;

for deg in range(1,50,1):
    i = 0
    j = 0

    while i < 2 * DegToStepInner:
        StepMotorOuter(CW)
        i+=1
        
    while j < DegToStepInner:
        StepMotorInner(CW)
        j+=1
    sleep(0.1)
    #data.append(myADC.readVoltage())
    #data2.append(myADC.readVoltage())
    #data3.append(myADC.readVoltage())
     
myDigital_IO.write(0b00000)         

# %%
angleInt = []
angleInt2 = []
angleInt3 = []

A = np.radians(45)
n = 1.52

def IntAngles(step):
    ang = np.radians(step / DegToStepInner)
    return  np.degrees(np.arcsin((np.sin(ang - A)/n)) + A)
    
for i in range(len(data)):
    angleInt.append(IntAngles(i+1))
    
for i in range(len(data2)):
    angleInt2.append(IntAngles(i+1))
    
for i in range(len(data3)):
    angleInt3.append(IntAngles(i+1))
    

# %%
print(data)
print(angleInt)

# %%
print(data2)
print(angleInt2)

# %%
print(data3)
print(angleInt3)

# %%
## Red laser
dataNorm = []

for i in data:
    dataNorm.append(i/max(data))

plt.plot(angleInt, dataNorm)
plt.grid()
plt.show()

# %%
## Green laser
data2Norm = []

for i in data2:
    data2Norm.append(i/max(data2))

plt.plot(angleInt2, data2Norm)
plt.grid()
plt.show()

# %%
## Blue laser
data3Norm = []

for i in data3:
    data3Norm.append(i/max(data3))

plt.plot(angleInt3, data3Norm)
plt.grid()
plt.show()