# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from ipywidgets import interact, interactive
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from IPython.display import display, Math, Latex

def balance_energia(Conversion=0, Expresion=1, Tref=0):
    # PARÁMETROS DEL SISTEMA
    Hf0A = -100
    Hf0B = -200
    Hf0C = -1000
    Hf0D = -500
    CpA = 5
    CpB = 2
    CpC = 1
    CpD = 1
    # Datos de input
    T1=0
    T2=50
    T3=100
    T4=75
    F1A = 10.
    F1D = 90.
    F2B = 15.
    F2D = 85.
    F3D = 50.
    # Balance de materia
    xi=F1A*Conversion
    F4A = F1A-xi
    F4B = F2B-xi
    F4C = 2*xi
    F4D = F1D+F2D-F3D
    
    # Balance de energia - Términos de calor sensible
    NcalCp=0
    NcalHf=0
    H1s = (F1A*CpA + F1D*CpD)*(T1-Tref)
    if H1s!=0:
        NcalCp=NcalCp+2
    H2s = (F2B*CpB + F2D*CpD)*(T2-Tref)
    if H2s!=0:
        NcalCp=NcalCp+2
    H3s = F3D*CpD*(T3-Tref)
    if H3s!=0:
        NcalCp=NcalCp+1
    H4s = (F4A*CpA + F4B*CpB + F4C*CpC + F4D*CpD)*(T4-Tref)
    if H4s!=0:
        if F4A==0 or F4C==0:
            NcalCp=NcalCp+3
        else:
            NcalCp=NcalCp+4
    Hr0 = xi*(2*Hf0C-Hf0A-Hf0B)
    Qreac = Hr0 + H3s + H4s - H1s - H2s
    if Expresion == 1:
        H1 = H1s
        H2 = H2s
        H3 = H3s
        H4 = H4s
        if Conversion != 0:
            txt0 = 'INCORRECT EXPRESSION FOR A BALANCE WITH CHEMICAL REACTION!!!'
        elif Tref != T4:
            txt0 = 'The chosen reference temperature leads to a higher number of calculations than the minimum necessary'
        else:
            txt0 = 'WELL-SELECTED EXPRESSION AND REFERENCE TEMPERATURE'
    else:
        H1 = F1A*Hf0A + F1D*Hf0D + H1s 
        H2 = F2B*Hf0B + F2D*Hf0D + H2s 
        H3 = F3D*Hf0D + H3s 
        H4 = F4A*Hf0A + F4B*Hf0B + F4C*Hf0C + F4D*Hf0D + H4s 
        if F4A==0 or F4C==0:
            NcalHf=8
        else:
            NcalHf=7
        if Tref != 25:
            txt0 = 'REFERENCE TEMPERATURE VALUE INCOMPATIBLE WITH THE USE OF STANDARD FORMATION HEATS!!!'
        elif Conversion != 0:
            txt0 = 'WELL-SELECTED EXPRESSION AND REFERENCE TEMPERATURE'
        else:
            txt0 = 'Balance without chemical reaction: The chosen expression leads to an exaggerated number of calculations'
    Q = H3+H4-H1-H2
    
    fig, ax1 = plt.subplots(figsize=(9,9))
    img = mpimg.imread('./images/pq_balance_energia.png')
    imgplot = plt.imshow(img)
    plt.axis('off')
    txt1 = "F$_{1A}$: "+str(int(F1A))+"\nF$_{1D}$: "+str(int(F1D))
    plt.text(-150,100,txt1, size=14)
    txt2 = "F$_{2B}$: "+str(int(F2B))+"\nF$_{2D}$: "+str(int(F2D))
    plt.text(-150,500,txt2, size=14)
    txt3 = "F$_{3D}$: "+str(int(F3D))
    plt.text(1750,100,txt3, size=14)
    txt4 = "F$_{4A}$: "+str(int(F4A))+"\nF$_{4B}$: "+str(int(F4B))+"\nF$_{4C}$: "+str(int(F4C))+"\nF$_{4D}$: "+str(int(F4D))
    plt.text(1750,500,txt4, size=14)
    txt5 = r'$\xi$: '+str(int(xi))
    plt.text(750,500,txt5, size=14)
    
    plt.text(0,750,txt0, size=18, color='red')
    txtH1 = 'Enthalpy of stream 1: H$_1$ = '+str(int(H1))+'\n'
    txtH2 = 'Enthalpy of stream 2: H$_2$ = '+str(int(H2))+'\n'
    txtH3 = 'Enthalpy of stream 3: H$_3$ = '+str(int(H3))+'\n'
    txtH4 = 'Enthalpy of stream 4: H$_4$ = '+str(int(H4))+'\n'
    txtQ = 'HEAT: Q = '+str(int(Q))+'\n'
    txtHr1 = r'$\xi \cdot \Delta$ H$_{Rx}^0$ = '+str(int(Hr0))+'\n'
    txtHr2 = 'Heat calculated using the heat of reaction method: '+str(int(Qreac))
    txtH=txtH1+txtH2+txtH3+txtH4+txtQ+txtHr1+txtHr2
    plt.text(0,1400,txtH, size=16)
    txtNcal1 = r'Number of Cp $\Delta$T terms calculated: '+str(int(NcalCp))+'\n'
    txtNcal2 = r'Number of $\Delta H_f$ terms included:' +str(int(NcalHf))
    txtNcal=txtNcal1+txtNcal2
    plt.text(0,1700,txtNcal, size=16)
    plt.show()
    
#w = interactive(balance_energia, {'manual': True}, Conversion=(0.,1.,0.5), Expresion=(1, 2), Tref=(0.,100.,25.))
w = interactive(balance_energia, {'manual': True}, Conversion=[0,0.5,1], Expresion=[1, 2], Tref=[('T0 = 25', 25), ('T1 = 0', 0), ('T2 = 50', 50), ('T3 = 100', 100), ('T4 = 75', 75)])

