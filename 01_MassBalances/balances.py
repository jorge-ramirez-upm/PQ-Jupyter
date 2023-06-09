from ipywidgets import interactive, fixed, FloatSlider, interactive
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import sympy as sym

def selectividad(F_Ain=5, F_Bin = 0, F_Cin = 0, F_Aout = 4, F_Bout = 0):
    
    # CONTROL INPUT CORRECTO
    if F_Aout > F_Ain:
        print("INPUT ERROR, INPUT REACTANT SHOULD BE GREATER THAN ITS OUTPUT")
        return
    
    if F_Bout < F_Bin:
        print("INPUT ERROR, INPUT B SHOULD BE GREATER THAN ITS OUTPUT")
        return
    
    if F_Aout-F_Ain + 0.5*(F_Bout-F_Bin) > 0:
        print("INPUT ERROR, B CANNOT GROW MORE THAN TWICE AS MUCH AS A DECREASES A")
        return
    
    # CALCULOS
    xi1 = (F_Bout - F_Bin)/2.
    deltaA = F_Ain - F_Aout
    xi2 = deltaA - xi1
    F_Cout = F_Cin + xi2
    
    # Conversion
    X = 100*deltaA/F_Ain
    # Selectividad
    S = 100*xi1/deltaA
    # Rendimiento
    R = 100*xi1/F_Ain
    
    fig, ax1 = plt.subplots(figsize=(7,7))
    img = mpimg.imread('./images/esquema_selectividad.png')
    imgplot = plt.imshow(img)
    plt.axis('off')
    txt1 = "A: "+str(round(F_Ain,2))+"\nB: "+str(round(F_Bin,2))+"\nC: "+str(round(F_Cin,2))
    plt.text(-200,200,txt1, size=20)
    txt2 = "A: "+str(round(F_Aout,2))+"\nB: "+str(round(F_Bout,2))+"\nC: "+str(round(F_Cout,2))
    plt.text(1050,200,txt2, size=20)
    
    txt00 = r'$\xi_1$: '+str(round(xi1,2))+'\n'+r'$\xi_2$: '+str(round(xi2,2))+'\n$X$ (%): '+str(round(X,2))+'\n$S$ (%): '+str(round(S,2))+'\n$\eta$ (%): '+str(round(R,2))
    plt.text(1500,250,txt00, size=20,bbox={'alpha': 0.5, 'pad': 10})
    plt.show()
    
w = interactive(selectividad, {'manual': True}, F_Ain=(0.1,10.), F_Bin=(0., 10.), F_Cin=(0., 10.), F_Aout=(0.,10.), F_Bout=(0.,20.))

def metanol(F6=6.25):
    # DATOS
    F1 = 100.        # Caudal de alimentación fresca
    Xpaso = 18.      # Conversión por paso en %
    x1_H2 = 67.3     # Fracción molar (%) de hidrógeno en la alimentación fresca
    x1_CO = 32.5     # Fracción molar (%) de monóxido de carbono en la alimentación fresca
   
    x1_CH4 = 100-x1_H2-x1_CO
    
    F6min = F1*(x1_CH4+x1_H2-2*x1_CO)/100.
    
    F6max = F1*(x1_CH4+x1_H2+x1_CO*(1-3*Xpaso/100.))/100.
    
    if F6 <= F6min or F6 >= F6max:
        print("PURGE FLOW RATE OUT OF ALLOWABLE VALUES")
    else:
    
        # Variables
        variables = sym.var('F1_H2 F1_CO F1_CH4 F2_H2 F2_CO F2_CH4 F3_H2 F3_CO F3_CH4 F3_CH3OH F4_CH3OH \
                F5_H2 F5_CO F5_CH4 F6_H2 F6_CO F6_CH4 F7_H2 F7_CO F7_CH4 xi')
    
        # Ecuaciones de balance
        balance_mezclador = [
            sym.Eq(F1_H2+F7_H2, F2_H2),
            sym.Eq(F1_CO+F7_CO, F2_CO),
            sym.Eq(F1_CH4+F7_CH4, F2_CH4)
        ]

        balance_reactor = [
            sym.Eq(F2_H2 - 2*xi, F3_H2),
            sym.Eq(F2_CO - xi, F3_CO),
            sym.Eq(F2_CH4, F3_CH4),
            sym.Eq(xi, F3_CH3OH)
        ]

        balance_separador = [
            sym.Eq(F3_H2, F5_H2),
            sym.Eq(F3_CO, F5_CO),
            sym.Eq(F3_CH4, F5_CH4),
            sym.Eq(F3_CH3OH, F4_CH3OH)
        ]

        balance_divisor = [
            sym.Eq(F5_H2, F6_H2 + F7_H2),
            sym.Eq(F5_CO, F6_CO + F7_CO),
            sym.Eq(F5_CH4, F6_CH4 + F7_CH4)
        ]

        espec = [
            sym.Eq(F1_H2, x1_H2*F1/100.),
            sym.Eq(F1_CO, x1_CO*F1/100.),
            sym.Eq(F1_CH4, x1_CH4*F1/100.),
            sym.Eq(F6, F6_H2 + F6_CO + F6_CH4),
            sym.Eq(F3_CO, (1-Xpaso/100.)*F2_CO),
            sym.Eq(F5_H2/F5_CO, F7_H2/F7_CO),
            sym.Eq(F5_CH4/F5_CO, F7_CH4/F7_CO)
        ]

        ecuaciones = balance_mezclador + balance_reactor + balance_separador + balance_divisor + espec

        soluc = sym.solve(ecuaciones)
        soln = soluc[0]
    
        x6_CH4 = 100.*soln[F6_CH4]/F6
        F2 = soln[F2_CO]+soln[F2_H2]+soln[F2_CH4]
        x2_CO = 100.*soln[F2_CO]/F2
        Xglobal = 100.*(soln[F1_CO]-soln[F6_CO])/soln[F1_CO]
        

        fig, ax1 = plt.subplots(figsize=(12,12))
        img = mpimg.imread('./images/esquema_recirculacion_1.png')
        imgplot = plt.imshow(img)
        plt.axis('off')
        txt1 = "H$_2$: "+str(round(soln[F1_H2],2))+"\nCO: "+str(round(soln[F1_CO],2))+"\nCH$_4$: "+str(round(soln[F1_CH4],2))
        plt.text(0,100,txt1, size=14)
        txt2 = "H$_2$: "+str(round(soln[F2_H2],2))+"\nCO: "+str(round(soln[F2_CO],2))+"\nCH$_4$: "+str(round(soln[F2_CH4],2))
        plt.text(350,100,txt2, size=14)
        txt3 = "H$_2$: "+str(round(soln[F3_H2],2))+"\nCO: "+str(round(soln[F3_CO],2))+"\nCH$_4$: "+str(round(soln[F3_CH4],2))+"\nCH$_3$OH: "+str(round(soln[F3_CH3OH],2))
        plt.text(1150,-20,txt3, size=14)
        txt4 = "CH$_3$OH: "+str(round(soln[F4_CH3OH],2))
        plt.text(2000,100,txt4, size=14)
        txt5 = "H$_2$: "+str(round(soln[F5_H2],2))+"\nCO: "+str(round(soln[F5_CO],2))+"\nCH$_4$: "+str(round(soln[F5_CH4],2))
        plt.text(1750,600,txt5, size=14)
        txt6 = "H$_2$: "+str(round(soln[F6_H2],2))+"\nCO: "+str(round(soln[F6_CO],2))+"\nCH$_4$: "+str(round(soln[F6_CH4],2))
        plt.text(2150,900,txt6, size=14)
        txt7 = "H$_2$: "+str(round(soln[F7_H2],2))+"\nCO: "+str(round(soln[F7_CO],2))+"\nCH$_4$: "+str(round(soln[F7_CH4],2))
        plt.text(350,800,txt7, size=14)
        
        txt00 = "OVERALL CONVERSION (%):                              "+str(round(Xglobal,2))+"\nReactor inlet flow rate (mol/h):            "+str(round(F2,2))+"\nCO content at reactor inlet (%): "+str(round(x2_CO,2))
        plt.text(2350,500,txt00, size=14,bbox={'alpha': 0.5, 'pad': 10})
        plt.show()

w2=interactive(metanol, {'manual': True}, F6=FloatSlider(min=2.60, max=82.40, value=6.25, description='Purge'))