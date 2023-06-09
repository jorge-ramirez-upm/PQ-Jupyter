import numpy as np
#%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib.image as image
from matplotlib.patches import Rectangle, Polygon
from scipy.optimize import fsolve
from ipywidgets import interact, interactive, fixed, interact_manual, widget, widgets, Layout, HBox, VBox
from IPython.display import Image
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as mpl_patches

# Create a dictionary with the value of Henry's constant for several gases in water at 298.15 K
# Data taken from (in mol/m3.Pa, transformed into atm^-1 by multiplying by the factor 1.83089
Henry = dict()
#Henry['Acetona'] = 1.83089*2.8e-1
Henry['O2'] =  1.0/4.3e4
Henry['H2'] = 1.0/7.1e4
Henry['CO2'] = 1.0/1.6e3
Henry['N2'] = 1.0/9.1e4
Henry['He'] = 1.0/1.5e5
Henry['Ne'] = 1.0/1.2e5
Henry['Ar'] = 1.0/4.0e4
Henry['CO'] = 1.0/5.8e4

# Solution enthalpy divided by R (gas constant) in K

HsolR = dict()
#HsolR['Acetona'] =  -5000
HsolR['O2'] =  -1700
HsolR['H2'] = -500
HsolR['CO2'] = -2400
HsolR['N2'] = -1300
HsolR['He'] = -230
HsolR['Ne'] = -490
HsolR['Ar'] = -1300
HsolR['CO'] = -1300

def HenryT(Gas='O2', T=298.15):
    Href=Henry[Gas]
    Hsol=HsolR[Gas]
    return Href*np.exp(-Hsol*(1/T-1/298.15))

def checksolubility(GasA='O2', GasB='N2', T=298.15):
    p = np.linspace(0,5,100)
    plt.figure(2,figsize=(6, 4), dpi= 100)
    ax = plt.subplot(121)
    xA = HenryT(GasA, T)*p
    index=xA<1
    plt.plot(p[index],xA[index],color='black',label=GasA)
    xB = HenryT(GasB, T)*p
    index=xB<1
    plt.plot(p[index],xB[index],color='red',label=GasB)
    plt.xlabel('$p_i$ [atm]')
    plt.ylabel('$x$ [-]')
    plt.legend()
    plt.minorticks_on()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1g'))
    #plt.grid(linewidth=1, which='both')

    ax = plt.subplot(122)
    Tv=np.linspace(150, 400, 100)
    plt.plot(Tv,HenryT(GasA, Tv),color='black',label=GasA)
    plt.plot(Tv,HenryT(GasB, Tv),color='red',label=GasB)
    plt.xlabel('T [K]')
    plt.ylabel('$H_i$ [$\mathrm{atm}^{-1}$]')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1g'))
    ax.yaxis.tick_right()
        
    plt.show()
    
widgetsolubility=interactive(checksolubility,
                             GasA=Henry.keys(),
                             GasB=Henry.keys(),
                             T=widgets.BoundedFloatText(value=298.15,min=150,max=373.15,step=5,description='T (K)')
)
controlssolubility = HBox(widgetsolubility.children[:-1], layout = Layout(flex_flow='row wrap'))
outputsolubility = widgetsolubility.children[-1]

def checksolubilityYX(GasA='O2', GasB='N2', T=298.15, P=1):
    plt.figure(2,figsize=(6, 4), dpi= 100)
    ax = plt.subplot(111)

    HA = HenryT(GasA, T)
    XmaxA = HA*P/(1-HA*P)/2
    X = np.linspace(0,XmaxA,100)
    YA = X/(HA*P*(1+X)-X)
    plt.plot(X,YA,color='black',label=GasA)

    HB = HenryT(GasB, T)
    XmaxB = HB*P/(1-HB*P)/2
    X = np.linspace(0,XmaxB,100)
    YB = X/(HB*P*(1+X)-X)
    plt.plot(X,YB,color='red',label=GasB)

    plt.xlabel('$X$ [-]')
    plt.ylabel('$Y$ [-]')
    plt.legend()
    plt.minorticks_on()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    #plt.grid(linewidth=1, which='both')
    
    plt.show()
    
widgetsolubilityYX=interactive(checksolubilityYX,
                             GasA=Henry.keys(),
                             GasB=Henry.keys(),
                             T=widgets.BoundedFloatText(value=298.15,min=160,max=373.15,step=5,description='T (K)'),
                             P=widgets.BoundedFloatText(value=1.0,min=0.2,max=10.0,step=0.1,description='P (atm)')
)
controlssolubilityYX = HBox(widgetsolubilityYX.children[:-1], layout = Layout(flex_flow='row wrap'))
outputsolubilityYX = widgetsolubilityYX.children[-1]

def checkYXmin(Y0=0.01, YL=0.1, X0=0.0):
    GasA='CO2'
    plt.figure(2,figsize=(6, 4), dpi= 100)
    ax = plt.subplot(111)

    T=298.15
    P=1
    HA = HenryT(GasA, T)
    Xs = YL*HA*P/(1-YL*HA*P+YL)

    X = np.linspace(0,1.25*Xs,100)
    YA = X/(HA*P*(1+X)-X)
    plt.plot(X,YA,color='black',label=GasA)

    Xs = YL*HA*P/(1-YL*HA*P+YL)
    plt.plot([X0, Xs], [Y0, YL], color='red', label="$(L'/G')_{min}$")
    
    plt.xlabel('$X$ [-]')
    plt.ylabel('$Y$ [-]')
    plt.legend()
    plt.minorticks_on()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    #plt.grid(linewidth=1, which='both')
    #plt.xlim([0, 1.25*Xs])
    plt.xlim([0, 1.2e-4])
    plt.ylim([0, np.max(YA)])
    ax = plt.axis()
    
    # TEXTS
    plt.text(X0,ax[2],'$X_0$',size=10, horizontalalignment='left', verticalalignment='top', color='g')
    plt.plot([X0,X0],[0,Y0],linestyle='--',color='g')
    plt.text(ax[0], Y0, '$Y_0$',size=10, horizontalalignment='right', verticalalignment='bottom', color='b')
    plt.plot([0,X0],[Y0,Y0],linestyle='--',color='b')

    plt.text(Xs,ax[2],'$X_L^*$',size=10, horizontalalignment='left', verticalalignment='top', color='g')
    plt.plot([Xs,Xs],[0,YL],linestyle='--',color='g')
    plt.text(ax[0], YL, '$Y_L$',size=10, horizontalalignment='right', verticalalignment='bottom', color='b')
    plt.plot([0,Xs],[YL,YL],linestyle='--',color='b')
    
    LpGpmin = (YL-Y0)/(Xs-X0)
    plt.text((X0+Xs)/2, (Y0+YL)/2, "$\\left( \\frac{L'}{G'} \\right) _{min}=%.3g$"%LpGpmin,size=10, horizontalalignment='right', verticalalignment='bottom', color='r')
    
    if X0/(HA*P*(1+X0)-X0) >= Y0:
        plt.text(ax[1]*0.9, Y0, 'Absorption is not possible',size=15, horizontalalignment='right', verticalalignment='bottom', color='b', bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))
    
    plt.show()
    
widgetYXmin=interactive(checkYXmin,
                        Y0=widgets.BoundedFloatText(value=0.01,min=0.001,max=0.05,step=0.001),
                        YL=widgets.BoundedFloatText(value=0.1,min=0.05,max=0.2,step=0.001),
                        X0=widgets.BoundedFloatText(value=0.0,min=0.0,max=0.00001,step=0.000001),
)
controlsYXmin = HBox(widgetYXmin.children[:-1], layout = Layout(flex_flow='row wrap'))
outputYXmin = widgetYXmin.children[-1]

def checkYXopera(Y0=0.01, YL=0.1, X0=0.0, F=1.5):
    GasA='CO2'
    plt.figure(2,figsize=(6, 4), dpi= 100)
    ax = plt.subplot(111)

    T=298.15
    P=1
    HA = HenryT(GasA, T)
    Xs = YL*HA*P/(1-YL*HA*P+YL)

    X = np.linspace(0,1.25*Xs,100)
    YA = X/(HA*P*(1+X)-X)
    plt.plot(X,YA,color='black',label=GasA)

    Xs = YL*HA*P/(1-YL*HA*P+YL)
    plt.plot([X0, Xs], [Y0, YL], color='red', label="$(L'/G')_{min}$")
    LpGpmin = (YL-Y0)/(Xs-X0)    
    LpGp=LpGpmin*F
    Xop=Xs/F
    plt.plot([X0, Xop], [Y0, YL], color='m', label="$(L'/G')$")
    
    plt.xlabel('$X$ [-]')
    plt.ylabel('$Y$ [-]')
    plt.legend()
    plt.minorticks_on()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    #plt.grid(linewidth=1, which='both')
    #plt.xlim([0, 1.25*Xs])
    plt.xlim([0, 1.2e-4])
    plt.ylim([0, np.max(YA)])
    ax = plt.axis()
        
    # TEXTS
    plt.text(X0,ax[2],'$X_0$',size=10, horizontalalignment='left', verticalalignment='bottom', color='g')
    plt.plot([X0,X0],[0,Y0],linestyle='--',color='g')
    plt.text(ax[0], Y0, '$Y_0$',size=10, horizontalalignment='right', verticalalignment='bottom', color='b')
    plt.plot([0,X0],[Y0,Y0],linestyle='--',color='b')

    plt.text(Xs,ax[2],'$X_L^*$',size=10, horizontalalignment='left', verticalalignment='bottom', color='g')
    plt.plot([Xs,Xs],[0,YL],linestyle='--',color='g')
    plt.text(ax[0], YL, '$Y_L$',size=10, horizontalalignment='right', verticalalignment='bottom', color='b')
    plt.plot([0,Xs],[YL,YL],linestyle='--',color='b')

    plt.text(Xop,ax[2],'$X_L$',size=10, horizontalalignment='left', verticalalignment='bottom', color='g')
    plt.plot([Xop,Xop],[0,YL],linestyle='--',color='g')    
    
    plt.text((X0+Xs)/2*1.3, (Y0+YL)/2, "$\\left( \\frac{L'}{G'} \\right) _{min}$=\n%.3g"%LpGpmin,size=10, horizontalalignment='center', verticalalignment='top', color='r')
    plt.text((X0+Xop)/2*0.7, (Y0+YL)/2, "$\\left( \\frac{L'}{G'} \\right)$=\n%.3g"%LpGp,size=10, horizontalalignment='center', verticalalignment='bottom', color='m')
    
    if X0/(HA*P*(1+X0)-X0) >= Y0:
        plt.text(ax[1]*0.9, Y0, 'Absorption is not possible',size=15, horizontalalignment='right', verticalalignment='bottom', color='b', bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))
    
    plt.show()
    
widgetYXopera=interactive(checkYXopera,
                        Y0=widgets.BoundedFloatText(value=0.01,min=0.001,max=0.05,step=0.001),
                        YL=widgets.BoundedFloatText(value=0.1,min=0.05,max=0.2,step=0.001),
                        X0=widgets.BoundedFloatText(value=0.0,min=0.0,max=0.00001,step=0.000001),
                        F=widgets.BoundedFloatText(value=1.5,min=1.1,max=5,step=0.1)
)
controlsYXopera = HBox(widgetYXopera.children[:-1], layout = Layout(flex_flow='row wrap'))
outputYXopera = widgetYXopera.children[-1]

def checkYXstep(Y0=0.01, YL=0.1, X0=0.0, F=1.5, Steps=1):
    GasA='CO2'
    plt.figure(2,figsize=(6, 4), dpi= 100)
    ax = plt.subplot(111)

    T=298.15
    P=1
    HA = HenryT(GasA, T)
    Xs = YL*HA*P/(1-YL*HA*P+YL)

    X = np.linspace(0,1.25*Xs,100)
    YA = X/(HA*P*(1+X)-X)
    plt.plot(X,YA,color='black',label=GasA)

    Xs = YL*HA*P/(1-YL*HA*P+YL)
    plt.plot([X0, Xs], [Y0, YL], color='red', linewidth=1, label="$(L'/G')_{min}$")
    LpGpmin = (YL-Y0)/(Xs-X0)    
    LpGp=LpGpmin*F
    XL=Xs/F
    plt.plot([X0, XL], [Y0, YL], color='m', label="$(L'/G')$")
    
    plt.xlabel('$X$ [-]')
    plt.ylabel('$Y$ [-]')
    #plt.legend()
    plt.minorticks_on()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    #plt.grid(linewidth=1, which='both')
    #plt.xlim([0, 1.25*Xs])
    plt.xlim([0, 1.2e-4])
    plt.ylim([0, np.max(YA)])
    ax = plt.axis()
        
    # TEXTS
    plt.text(X0,ax[2],'$X_0$',size=10, horizontalalignment='left', verticalalignment='bottom', color='g')
    plt.plot([X0,X0],[0,Y0],linestyle='--',color='g')
    plt.text(ax[0], Y0, '$Y_0$',size=10, horizontalalignment='right', verticalalignment='bottom', color='b')
    plt.plot([0,X0],[Y0,Y0],linestyle='--',color='b')

    plt.text(Xs,ax[2],'$X_L^*$',size=10, horizontalalignment='left', verticalalignment='bottom', color='g')
    plt.plot([Xs,Xs],[0,YL],linestyle='--',color='g')
    plt.text(ax[0], YL, '$Y_L$',size=10, horizontalalignment='right', verticalalignment='bottom', color='b')
    plt.plot([0,Xs],[YL,YL],linestyle='--',color='b')

    plt.text(XL,ax[2],'$X_L$',size=10, horizontalalignment='left', verticalalignment='bottom', color='g')
    plt.plot([XL,XL],[0,YL],linestyle='--',color='g')    
        
    Possible=True    
    if X0/(HA*P*(1+X0)-X0) >= Y0:
        plt.text(ax[1]*0.9, Y0, 'Absorption is not possible',size=15, horizontalalignment='right', verticalalignment='bottom', color='b', bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))
        Possible=False
    
    # START solution
    XP = Xs
    YP = YL
    nTray = 0

    while YP > Y0 and nTray < Steps and Possible:
        nTray += 1
        XQ = XP   
        
        XP = X0 + (YP-Y0)*(XL-X0)/(YL-Y0)
        plt.plot([XQ,XP],[YP,YP],'y')
        plt.text(XP*0.9,YP,nTray)

        YQ = YP
        YP = XP/(HA*P*(1+XP)-XP) # Equilibrium curve
        plt.plot([XP,XP],[YQ,YP],'y')
        
    if YP <= Y0:
        # The following is to use the legend to place the text in a smart way
        
        # create a list with two empty handles (or more if needed)
        handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                         lw=0, alpha=0)] * 3

        # create the corresponding number of labels (= the text you want to display)
        labels = []
        labels.append("%d eq. stages"%nTray)
        labels.append("$(L'/G')_{min}$=%.3g"%LpGpmin)
        labels.append("$(L'/G')$=%.3g"%LpGp)

        # create the legend, supressing the blank space of the empty line symbol and the
        # padding between symbol and label by setting handlelenght and handletextpad
        legend_properties = {'weight':'bold'}
        plt.legend(handles, labels, loc='best',
                  fancybox=True, framealpha=0.7, 
                  handlelength=0, handletextpad=0, prop=legend_properties)        
        
    
    plt.show()
    
widgetYXstep=interactive(checkYXstep,
                        Y0=widgets.BoundedFloatText(value=0.01,min=0.001,max=0.05,step=0.001),
                        YL=widgets.BoundedFloatText(value=0.1,min=0.05,max=0.2,step=0.001),
                        X0=widgets.BoundedFloatText(value=0.0,min=0.0,max=0.00001,step=0.000001),
                        F=widgets.BoundedFloatText(value=1.5,min=1.1,max=5,step=0.1),
                        Steps=widgets.BoundedIntText(value=1,min=1,max=50,step=1)
)
controlsYXstep = HBox(widgetYXstep.children[:-1], layout = Layout(flex_flow='row wrap'))
outputYXstep = widgetYXstep.children[-1]

def checkYXfull(Gas='CO2', Y0=0.01, YL=0.1, X0=0.0, F=1.5, T=298.15, P=1):
    plt.figure(2,figsize=(6, 4), dpi= 100)
    ax = plt.subplot(111)

    HA = HenryT(Gas, T)
    Xs = YL*HA*P/(1-YL*HA*P+YL)

    X = np.linspace(0,1.25*Xs,100)
    YA = X/(HA*P*(1+X)-X)
    plt.plot(X,YA,color='black',label=Gas)

    Xs = YL*HA*P/(1-YL*HA*P+YL)
    plt.plot([X0, Xs], [Y0, YL], color='red', linewidth=1, label="$(L'/G')_{min}$")
    LpGpmin = (YL-Y0)/(Xs-X0)    
    LpGp=LpGpmin*F
    XL=Xs/F
    plt.plot([X0, XL], [Y0, YL], color='m', label="$(L'/G')$")
    
    plt.xlabel('$X$ [-]')
    plt.ylabel('$Y$ [-]')
    #plt.legend()
    plt.minorticks_on()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    #plt.grid(linewidth=1, which='both')
    plt.xlim([0, np.max(X)])
    plt.ylim([0, np.max(YA)])
    ax = plt.axis()
        
    # TEXTS
    plt.text(X0,ax[2],'$X_0$',size=10, horizontalalignment='left', verticalalignment='bottom', color='g')
    plt.plot([X0,X0],[0,Y0],linestyle='--',color='g')
    plt.text(ax[0], Y0, '$Y_0$',size=10, horizontalalignment='right', verticalalignment='bottom', color='b')
    plt.plot([0,X0],[Y0,Y0],linestyle='--',color='b')

    plt.text(Xs,ax[2],'$X_L^*$',size=10, horizontalalignment='left', verticalalignment='bottom', color='g')
    plt.plot([Xs,Xs],[0,YL],linestyle='--',color='g')
    plt.text(ax[0], YL, '$Y_L$',size=10, horizontalalignment='right', verticalalignment='bottom', color='b')
    plt.plot([0,Xs],[YL,YL],linestyle='--',color='b')

    plt.text(XL,ax[2],'$X_L$',size=10, horizontalalignment='left', verticalalignment='bottom', color='g')
    plt.plot([XL,XL],[0,YL],linestyle='--',color='g')    
        
    Possible=True    
    if X0/(HA*P*(1+X0)-X0) >= Y0:
        plt.text(ax[1]*0.9, Y0, 'Absorption is not possible',size=15, horizontalalignment='right', verticalalignment='bottom', color='b', bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))
        Possible=False
    
    # START solution
    XP = Xs
    YP = YL
    nTray = 0

    while YP > Y0 and Possible:
        nTray += 1
        XQ = XP   
        
        XP = X0 + (YP-Y0)*(XL-X0)/(YL-Y0)
        plt.plot([XQ,XP],[YP,YP],'y')
        plt.text(XP*0.9,YP,nTray)

        YQ = YP
        YP = XP/(HA*P*(1+XP)-XP) # Equilibrium curve
        plt.plot([XP,XP],[YQ,YP],'y')
        
    if YP <= Y0:
        # The following is to use the legend to place the text in a smart way
        
        # create a list with two empty handles (or more if needed)
        handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                         lw=0, alpha=0)] * 3

        # create the corresponding number of labels (= the text you want to display)
        labels = []
        labels.append("%d eq. stages"%nTray)
        labels.append("$(L'/G')_{min}$=%.3g"%LpGpmin)
        labels.append("$(L'/G')$=%.3g"%LpGp)

        # create the legend, supressing the blank space of the empty line symbol and the
        # padding between symbol and label by setting handlelenght and handletextpad
        legend_properties = {'weight':'bold'}
        plt.legend(handles, labels, loc='best',
                  fancybox=True, framealpha=0.7, 
                  handlelength=0, handletextpad=0, prop=legend_properties)        
        
    
    plt.show()
    
widgetYXfull=interactive(checkYXfull,
                        Gas=Henry.keys(),
                        Y0=widgets.BoundedFloatText(value=0.01,min=0.001,max=0.05,step=0.001),
                        YL=widgets.BoundedFloatText(value=0.1,min=0.05,max=0.2,step=0.001),
                        X0=widgets.BoundedFloatText(value=0.0,min=0.0,max=0.00001,step=0.00000001),
                        F=widgets.BoundedFloatText(value=1.5,min=1.1,max=5,step=0.1),
                        T=widgets.BoundedFloatText(value=298.15,min=273.15,max=373.15,step=5,description='T (K)'),
                        P=widgets.BoundedFloatText(value=1.0,min=0.2,max=10.0,step=0.1,description='P (atm)') 
)
controlsYXfull = HBox(widgetYXfull.children[:-1], layout = Layout(flex_flow='row wrap'))
outputYXfull = widgetYXfull.children[-1]
