import scipy.optimize as optimize
import matplotlib.pyplot as plt
import math
import numpy as np
import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk
import time
start = time.time()

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# File save dialog
def file_save(text):
    f = Tk.filedialog.asksaveasfile(mode='w', defaultextension=".fit")
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    f.write(text)
    f.close()

# File load dialog
def file_load():
    file_path = Tk.filedialog.askopenfilename()
    s = []
    I = []
    Err = []
    with open(file_path) as f:
        plotlist = f.readlines()
    for str in plotlist:
        if not str.startswith(' '): continue
        ss = str.split()
        if ss.__len__() != 3: continue
        if is_number(ss[0]):
            s.append(float(ss[0]))
        if is_number(ss[1]):
            I.append(float(ss[1]))
        if is_number(ss[2]):
            Err.append(float(ss[2]))
    return s, I, Err

def SAXSplot(s, I, Err):

    plt.figure()
    plt.subplot()
    # TODO: errors -> fit calculation (I/sigma)
    #plt.errorbar(s, I, yerr = Err, fmt=".-", color = 'blue')
    plt.scatter(s,I,s=80, facecolors='none', edgecolors='b')
    plt.yscale("log")
    plt.xlabel("s, $\mathregular{nm^{-1}}$")
    plt.ylabel("Intensity, a.u.")

def SAXSplotfit(q, I, Q1):
    plt.subplot()
    Err1,Err2 = [], []
    s1,s2 = [], []
    i1,i2 = [], []
    for sca,inten in zip(q,I):
        if sca <= Q1:
            s1.append(sca)
            i1.append(inten)
            Err1.append(math.sqrt(math.fabs(inten)))

        if sca >= Q1:
            s2.append(sca)
            i2.append(inten)
            Err2.append(math.sqrt(math.fabs(inten)))
    plt.plot(s1, i1, linewidth = 2, color = 'yellow')
    plt.plot(s2, i2, linewidth = 2, color = 'red')

#function to fit 6 params
def Guinier_Porod(q,G,s,d,Rg,b, Q1):
    result = []
    G*=1000
    b*=1000
    if (d - s)  < 0:
        return 0
    for qq in q:
        if qq <= Q1: result.append(b+((G/math.pow(qq,s))*np.exp(-(math.pow(qq*Rg,2))/(3-s))))
        if qq > Q1: result.append(b+((G/math.pow(qq,d))*math.pow(Q1,d-s)*np.exp(-(math.pow(Q1*Rg,2))/(3-s))))
    if (len(result) < len (q)) : return np.zeros(len(q))
    return np.array(result)

#function to fit 5 params
def Guinier_Porod5param(q,G,s,d,Rg,b):
    result = []
    G*=1000
    b*=1000
    if (d - s)  < 0:
        return 0
    Q1 = (1 / Rg) * math.pow((0.5 * (d - s) * (3 - s)), 0.5)
    for qq in q:
        if qq <= Q1: result.append(b+((G/math.pow(qq,s))*np.exp(-(math.pow(qq*Rg,2))/(3-s))))
        if qq > Q1: result.append(b+((G/math.pow(qq,d))*math.pow(Q1,d-s)*np.exp(-(math.pow(Q1*Rg,2))/(3-s))))
    if (len(result) < len (q)) : return np.zeros(len(q))
    return np.array(result)

#///////////////////////START PROGRAM////////////////////////////////////////////////////////////

# open and read pdb file
root = Tk.Tk()
root.withdraw()
q, I, Err = file_load()

# Create x for fit
trialQ = np.linspace(q[0], q[-1], 1000)
# Fit
popt, pcov = optimize.curve_fit(Guinier_Porod5param, q, I, bounds = ([1,0,1,0.4,0], [600,3,4,10.0,1]),maxfev=100000)
G,s,d,Rg,background = popt
Q1 = (1 / Rg) * math.pow((0.5 * (d - s) * (3 - s)), 0.5)
str = "G = " + str(G*1000) + " s = " + str(s) + " d = " + str(d) + " Rg = " + str(Rg) +" background = " +str(background*1000) +\
      "Q1 = " + str(Q1)
print(str)
yEXP = Guinier_Porod5param(trialQ, *popt)
SAXSplot(q,I,Err)
SAXSplotfit(trialQ, yEXP, Q1)
plt.title("Generalized Guinier-Porod fit \n" + str )
plt.grid()
plt.grid(b=True, which='minor', color='black', linestyle='--')
plt.show()
end = time.time()
print (" Time consumed : " + str (end - start) + " sec")

