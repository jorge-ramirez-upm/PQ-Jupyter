import numpy as np
#%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib.image as image
from matplotlib.patches import Rectangle, Polygon
from scipy.optimize import fsolve
from ipywidgets import interact, interactive, fixed, interact_manual, widget, widgets, Layout, HBox, VBox
from IPython.display import Image

# Create a dictionary with the Antoine Equation for each of the components
# We use a dictionary to define the Antoine functions for a series of substances. 
# Data from http://ddbonline.ddbst.com/AntoineCalculation/AntoineCalculationCGI.exe
Psat = dict()
Psat['Acetone'] = lambda T: 10**(7.1327 - 1219.97/(T + 230.653))
Psat['Acetonitrile'] = lambda T: 10**(7.33986 - 1482.29/(T + 250.523))
Psat['Acetic acid'] = lambda T: 10**(7.2996 - 1479.02/(T + 216.82))
Psat['Water'] = lambda T: 10**(8.07131 - 1730.63/(T + 233.426))
Psat['Ethanol'] = lambda T: 10**( 8.20417 - 1642.89/(T + 230.3))
#Psat['Etilenglicol'] = lambda T: 10**( 8.7945 - 2615.4/(T + 244.91))
Psat['Phenol'] = lambda T: 10**( 7.1345 - 1516.07/(T + 174.57))
Psat['Isopropyl-alcohol'] = lambda T: 10**( 8.1182 - 1580.92/(T + 219.62))

Molecule = dict()
Molecule['Acetone'] = 'img/acetone.png'
Molecule['Acetonitrile'] = 'img/acetonitrile.png'
Molecule['Acetic acid'] = 'img/aceticacid.png'
Molecule['Water'] = 'img/water.png'
Molecule['Ethanol'] = 'img/ethanol.png'
#Molecule['Etilenglicol'] = 'img/ethylenglycol.png'
Molecule['Phenol'] = 'img/phenol.png'
Molecule['Isopropyl-alcohol'] = 'img/isopropylalcohol.png'

Hvap = dict()
Hvap['Acetone'] = 31300
Hvap['Acetonitrile'] = 33225
Hvap['Acetic acid'] = 23700
Hvap['Water'] = 40660
Hvap['Ethanol'] = 38600
#Hvap['Etilenglicol'] = 53200
Hvap['Phenol'] = 58800
Hvap['Isopropyl-alcohol'] = 44000

cPLiq = dict()
cPLiq['Acetone'] = 125.5
cPLiq['Acetonitrile'] = 91.7
cPLiq['Acetic acid'] = 123.1
cPLiq['Water'] = 75.327
cPLiq['Ethanol'] = 112.4
#cPLiq['Etilenglicol'] = 37.8
cPLiq['Phenol'] = 138.0
cPLiq['Isopropyl-alcohol'] = 160.8 

cPVap = dict()
cPVap['Acetone'] = 75
cPVap['Acetonitrile'] = 80 
cPVap['Acetic acid'] = 63.4
cPVap['Water'] = 37.47
cPVap['Ethanol'] = 82.0
#cPVap['Etilenglicol'] = 23.0
cPVap['Phenol'] = 145.52
cPVap['Isopropyl-alcohol'] = 92.4

def checkvolatility(FluidoA='Acetone', FluidoB='Water'):
    T = np.linspace(0,150,100)
    plt.figure(2,figsize=(6, 4), dpi= 100)    
    plt.plot(T,Psat[FluidoA](T)/760,color='black',label=FluidoA)
    plt.plot(T,Psat[FluidoB](T)/760,color='red',label=FluidoB)
    plt.xlabel('T $\circ$C')
    plt.ylabel('$p^0$ [atm]')
    plt.legend()
    plt.minorticks_on()
    #plt.grid(linewidth=1, which='both')
    
    xrng=plt.xlim()
    yrng=plt.ylim()
    im =image.imread(Molecule[FluidoA])
    scale=.3 #the image takes this fraction of the graph 
    scaley = scale*im.shape[0]/im.shape[1]
    plt.imshow(im,aspect='auto',extent=(xrng[0],xrng[0] + scale*(xrng[1]-xrng[0]), yrng[1] - (scale+scaley)*(yrng[1]-yrng[0]), yrng[1] - scale*(yrng[1]-yrng[0])), zorder=-1)
    plt.xlim(xrng)
    plt.ylim(yrng)
    xrng=plt.xlim()
    yrng=plt.ylim()
    im2 =image.imread(Molecule[FluidoB])
    scale=.3 #the image takes this fraction of the graph
    scaley = scale*im2.shape[0]/im2.shape[1]
    plt.imshow(im2,aspect='auto',extent=(xrng[0] + scale*(xrng[1]-xrng[0]), xrng[0] + 2*scale*(xrng[1]-xrng[0]), yrng[1] - (scale+scaley)*(yrng[1]-yrng[0]), yrng[1] - scale*(yrng[1]-yrng[0])), zorder=-1)
    plt.xlim(xrng)
    plt.ylim(yrng)
    
    plt.show()
    
widgetvolatility=interactive(checkvolatility,
                             FluidoA=Psat.keys(),
                             FluidoB=Psat.keys()                   
)
controlsvolatility = HBox(widgetvolatility.children[:-1], layout = Layout(flex_flow='row wrap'))
outputvolatility = widgetvolatility.children[-1]

Tsat = dict()
for s in Psat.keys():
    Tsat[s] = lambda P, s=s: fsolve(lambda T: Psat[s](T)-P,50)[0]
        
def getTeb(Fluido='Water', P=1):
    print("Boiling point of %s : %.3f °C"%(Fluido,Tsat[Fluido](P*760)))
    from IPython.display import Image
    u=Image(filename=Molecule[Fluido], width=150, height=200)
    display(u)
    
widgetTeb=interactive(getTeb, Fluido=Tsat.keys(), P=widgets.BoundedFloatText(value=1,min=0.1,max=5.0,step=0.1,description='p (atm)'))
controlsTeb = HBox(widgetTeb.children[:-1], layout = Layout(flex_flow='row wrap'))
outputTeb = widgetTeb.children[-1]

def checkTeb(FluidoA='Acetone', FluidoB='Water'):
    P = np.linspace(100,5000,100)
    TA=np.zeros(len(P))
    TB=np.zeros(len(P))
    for i in range(len(P)):
        TA[i]=Tsat[FluidoA](P[i])
        TB[i]=Tsat[FluidoB](P[i])
    plt.figure(2,figsize=(6, 4), dpi= 100)    
    plt.plot(P/760,TA,color='black',label=FluidoA)
    plt.plot(P/760,TB,color='red',label=FluidoB)
    plt.xlabel('$p_{ext}$ [atm]')
    plt.ylabel('T [$\circ$C]')
    plt.legend()
    plt.minorticks_on()
    plt.grid(linewidth=1, which='both')
    plt.show()
    
widgetTeb2=interactive(checkTeb, FluidoA=Tsat.keys(), FluidoB=Tsat.keys())
controlsTeb2 = HBox(widgetTeb2.children[:-1], layout = Layout(flex_flow='row wrap'))
outputTeb2 = widgetTeb2.children[-1]

def checkRaoult(FluidoA='Acetone', FluidoB='Water',T=100):
    plt.figure(2,figsize=(6, 4), dpi= 100)    
    PA=Psat[FluidoA](T)
    PB=Psat[FluidoB](T)
    x=np.linspace(0,1,100)
    plt.plot(x,PA*x/760,color='blue',label="$p^0_{"+FluidoA+"}$")
    plt.plot(x,PB*(1-x)/760,color='red',label="$p^0_{"+FluidoB+"}$")
    plt.plot(x,PB*(1-x)/760+PA*x/760,color='black',label="$p^0_{TOTAL}$")
    plt.xlabel("$x_{"+FluidoA+"}$")
    plt.ylabel('$p$ [atm]')
    plt.legend()
    plt.minorticks_on()
    plt.grid(linewidth=1, which='both')
    plt.xlim(0,1)
    plt.ylim(0,1.1*max(PA,PB)/760)
    plt.show()
    
widgetRaoult=interactive(checkRaoult, FluidoA=Psat.keys(),FluidoB=Psat.keys(),T=widgets.BoundedFloatText(value=100,min=0,max=200,step=5,description='T (°C)'))
controlsRaoult = HBox(widgetRaoult.children[:-1], layout = Layout(flex_flow='row wrap'))
outputRaoult = widgetRaoult.children[-1]

def checkTx(FluidoA='Acetone', FluidoB='Water',Pext=1.0):
    if FluidoA==FluidoB:
        print("Choose two different fluids")
        return
    P=Pext*760
    x = lambda T: (P-Psat[FluidoB](T))/(Psat[FluidoA](T)-Psat[FluidoB](T))
    y = lambda T: x(T)*Psat[FluidoA](T)/P
    plt.figure(2,figsize=(6, 4), dpi= 100)
    T = np.linspace(Tsat[FluidoA](P),Tsat[FluidoB](P))
    plt.plot([x(T) for T in T],T,color='black')
    plt.plot([y(T) for T in T],T,color='black')
    plt.xlabel("$x_{"+FluidoA+"}$")
    plt.ylabel('Temperature $^\circ$C')
    plt.title('Txy diagram for {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
    Tmid=0.5*(Tsat[FluidoA](P)+Tsat[FluidoB](P))
    plt.text(x(Tmid),Tmid-5,'L',size=15)
    plt.text(y(Tmid),Tmid+5,'V',size=15)
    plt.minorticks_on()
    plt.grid(linewidth=1, which='both')
    plt.xlim(0,1)
    plt.show()
    
widgetTx=interactive(checkTx,FluidoA=Psat.keys(),FluidoB=Psat.keys(),Pext=widgets.BoundedFloatText(value=1.0,min=0.1,max=5.0,step=0.1,description='p (atm)'))
controlsTx = HBox(widgetTx.children[:-1], layout = Layout(flex_flow='row wrap'))
outputTx = widgetTx.children[-1]

def checkFlash(FluidoA='Acetone', FluidoB='Water',Pext=1.0,xF=0.5,TF=100.0,F=100.0):
    if FluidoA==FluidoB:
        print("Choose two different fluids")
        return
    P=Pext*760
    x = lambda T: (P-Psat[FluidoB](T))/(Psat[FluidoA](T)-Psat[FluidoB](T))
    y = lambda T: x(T)*Psat[FluidoA](T)/P
    plt.figure(2,figsize=(6, 4), dpi= 100)
    T = np.linspace(Tsat[FluidoA](P),Tsat[FluidoB](P))
    plt.plot([x(T) for T in T],T,color='black')
    plt.plot([y(T) for T in T],T,color='black')
    plt.xlabel("$x_{"+FluidoA+"}$")
    plt.ylabel('Temperature $^\circ$C')
    plt.title('Txy diagram for {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
    plt.minorticks_on()
    plt.grid(linewidth=1, which='both')
    plt.xlim(0,1)

    Tdew = fsolve(lambda T: y(T)-xF, 138)
    Tbub = fsolve(lambda T: x(T)-xF, 0.01)
    
    ax = plt.axis()
    plt.plot(xF,TF,'kx',ms=8)
    plt.plot([xF,xF,0],[ax[2],TF,TF],'b--')
    plt.text(xF,ax[2],'$x_F$',size=20)
    plt.text(0.01,TF+1,'$T_F$',size=15)
    
    plt.plot(xF,Tdew,'kD',ms=8)
    plt.plot(xF,Tbub,'kD',ms=8)
    
    if (TF>=Tbub and TF<=Tdew):
        plt.plot(y(TF),TF,'go',ms=7)
        plt.plot([xF,y(TF),y(TF)],[TF,TF,ax[2]],'g--')
        plt.text(y(TF),ax[2],'$y=${:.2f}'.format(y(TF)),size=15)
        plt.plot(x(TF),TF,'ro',ms=7)
        plt.plot([xF,x(TF),x(TF)],[TF,TF,ax[2]],'r--')
        plt.text(x(TF),ax[2],'$x=${:.2f}'.format(x(TF)),size=15,horizontalalignment='right')
        L=F*(y(TF)-xF)/(y(TF)-x(TF))
        V=F-L
        print("L = %.3f mol/s -- V = %.3f mol/s -- Tbubble = %.1f °C -- Tdew = %.f °C"%(L, V, Tbub[0], Tdew[0]))
    plt.show()
        
widgetFlash=interactive(checkFlash,
                        FluidoA=Psat.keys(),
                        FluidoB=Psat.keys(),
                        Pext=widgets.BoundedFloatText(value=1.0,min=0.1,max=5.0,step=0.1,description='p (atm)'),
                        xF=widgets.BoundedFloatText(value=0.5,min=0.05,max=0.95,step=0.01,description='$x_F$'),
                        TF=widgets.BoundedFloatText(value=100,min=0,max=200,step=5,description='$T_F$ (°C)'),
                        F=widgets.BoundedFloatText(value=10.0,min=1.0,max=100.0,step=1.0,description='F (mol/s)')
                       )
controlsFlash = HBox(widgetFlash.children[:-1], layout = Layout(flex_flow='row wrap'))
outputFlash = widgetFlash.children[-1]

def checkyx(FluidoA='Acetone', FluidoB='Water',Pext=1.0):
    if FluidoA==FluidoB:
        print("Choose two different fluids")
        return
    P=Pext*760
    x = lambda T: (P-Psat[FluidoB](T))/(Psat[FluidoA](T)-Psat[FluidoB](T))
    y = lambda T: x(T)*Psat[FluidoA](T)/P
    T = np.linspace(Tsat[FluidoA](P),Tsat[FluidoB](P))
    plt.figure(figsize=(5,5),dpi=100)
    plt.plot([x(T) for T in T],[y(T) for T in T], color='black')
    plt.plot([0,1],[0,1],color='black',linestyle='--')
    plt.axis('equal')

    plt.title('x-y diagram for {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
    plt.xlabel("$x_{"+FluidoA+"}$")
    plt.ylabel("$y_{"+FluidoA+"}$")

    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.minorticks_on()
    plt.grid(linewidth=1, which='both')

widgetyx=interactive(checkyx,
                     FluidoA=Psat.keys(),
                     FluidoB=Psat.keys(),
                     Pext=widgets.BoundedFloatText(value=1.0,min=0.1,max=5.0,step=0.1,description='p (atm)')
                     )
controlsyx = HBox(widgetyx.children[:-1], layout = Layout(flex_flow='row wrap'))
outputyx = widgetyx.children[-1]

##### Recta de alimentación
def checkFeed(TF=100.0,xF=0.5):
    FluidoA='Isopropyl-alcohol'
    FluidoB='Acetic acid'
    Pext=2.0
    
    P=Pext*760
    
    # Diatrama T-xy
    x = lambda T: (P-Psat[FluidoB](T))/(Psat[FluidoA](T)-Psat[FluidoB](T))
    y = lambda T: x(T)*Psat[FluidoA](T)/P
    plt.figure(2,figsize=(8, 4), dpi= 100)
    plt.subplot(121)
    T = np.linspace(Tsat[FluidoA](P),Tsat[FluidoB](P))
    plt.plot([x(T) for T in T],T,color='black')
    plt.plot([y(T) for T in T],T,color='black')
    plt.xlabel("$x_{"+FluidoA+"}$")
    plt.ylabel('Temperature $^\circ$C')
    plt.minorticks_on()
    plt.grid(linewidth=1, which='both')
    plt.xlim(0,1)
    plt.ylim(100,150)

    Tdew = fsolve(lambda T: y(T)-xF, 138)
    Tbub = fsolve(lambda T: x(T)-xF, 0.01)
    
    ax = plt.axis()
    plt.plot(xF,TF,'kx',ms=8)
    plt.plot([xF,xF,0],[ax[2],TF,TF],'b--')
    plt.text(xF+0.02,ax[2],'$x_F$',horizontalalignment='left',verticalalignment='bottom')
    plt.text(0.01,TF+1,'$T_F$')

    if (TF>=Tbub and TF<=Tdew):
        plt.plot(y(TF),TF,'go',ms=7)
        plt.plot([xF,y(TF),y(TF)],[TF,TF,ax[2]],'g--')
        plt.text(y(TF)+0.02,ax[2],'$y$',color='g',horizontalalignment='left',
                 verticalalignment='bottom')
        plt.plot(x(TF),TF,'ro',ms=7)
        plt.plot([xF,x(TF),x(TF)],[TF,TF,ax[2]],'r--')
        plt.text(x(TF)-0.02,ax[2],'$x$',horizontalalignment='right',
                 verticalalignment='bottom',color='r')
    
    # Diagrama y-x
    x = lambda T: (P-Psat[FluidoB](T))/(Psat[FluidoA](T)-Psat[FluidoB](T))
    y = lambda T: x(T)*Psat[FluidoA](T)/P
    T = np.linspace(Tsat[FluidoA](P),Tsat[FluidoB](P),100)
    #plt.figure(figsize=(5,5),dpi=100)
    plt.subplot(122)
    xT=x(T)
    yT=y(T)
    p=xT.argsort()
    xT=xT[p]
    yT=yT[p]
    plt.plot(xT,yT, color='black')
    plt.plot([0,1],[0,1],color='black',linestyle='--')
    plt.axis('equal')

    # RECTA Q
    Tdew = fsolve(lambda T: y(T)-xF, 138)[0]
    Tbub = fsolve(lambda T: x(T)-xF, 0.01)[0]
    
    Hvapmix=xF*Hvap[FluidoA]+(1-xF)*Hvap[FluidoB]
    cPVapmix=xF*cPVap[FluidoA]+(1-xF)*cPVap[FluidoB]
    cPLiqmix=xF*cPLiq[FluidoA]+(1-xF)*cPLiq[FluidoB]
    
    if (TF>Tdew):
        q = -(TF-Tdew)*(cPVapmix)/Hvapmix
    elif (TF>Tbub):
        q = (Tdew-TF)/(Tdew-Tbub)
    else:
        q = (Hvapmix+(Tbub-TF)*cPLiqmix)/Hvapmix
    
    yq=-xT*q/(1-q)+xF/(1-q)

    #Intersection occurs between points itemindex and itemindex+1
    itemindex = np.argwhere(np.diff(np.sign(yq - yT))).flatten()
    a1=np.array([xT[itemindex][0],yq[itemindex][0]])
    a2=np.array([xT[itemindex+1][0],yq[itemindex+1][0]])
    b1=np.array([xT[itemindex][0],yT[itemindex][0]])
    b2=np.array([xT[itemindex+1][0],yT[itemindex+1][0]])
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = np.empty_like(da)
    dap[0] = -da[1]
    dap[1] = da[0]
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    intersec=(num / denom.astype(float))*db + b1    
    plt.plot([xF,intersec[0]],[xF,intersec[1]],color='red')
    plt.plot(intersec[0],intersec[1],'co',ms=5)
    if (TF>=Tbub and TF<=Tdew):
        plt.text(intersec[0]-0.02,intersec[1]+0.02,'$(x,y)$',
                 horizontalalignment='right',verticalalignment='bottom')    
    
    plt.xlabel("$x_{"+FluidoA+"}$")
    plt.ylabel("$y_{"+FluidoA+"}$")

    plt.plot(xF,xF,'go',ms=5)

    
    plt.plot([xF,xF,0],[0,xF,xF],'b--')
    plt.text(xF+0.01,0.02,'$x_F$',horizontalalignment='left',verticalalignment='bottom')    
    plt.text(0.01, xF+0.01, '$x_F$',horizontalalignment='left',verticalalignment='bottom')    
    
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.minorticks_on()
    plt.grid(linewidth=1, which='both')

    plt.suptitle('T-xy and y-x diagrams for {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
    plt.subplots_adjust(wspace=0.2)
    plt.show()
    
widgetFeed=interactive(checkFeed,
                         TF=widgets.BoundedFloatText(value=100,min=100,max=150,step=1,description='$T_F$ (°C)'), 
                         xF=widgets.BoundedFloatText(value=0.5,min=0.1,max=0.9,step=0.01,description='$x_F$'))
controlsFeed = HBox(widgetFeed.children[:-1], layout = Layout(flex_flow='row wrap'))
outputFeed = widgetFeed.children[-1]

##### Diagrama y-x con recta de alimentación y recta de operación de la zona de rectificación

def checkFeedRectif(TF=100.0,xF=0.5,xD=0.9,Rext=3.0):
    FluidoA='Isopropyl-alcohol'
    FluidoB='Acetic acid'
    Pext=2.0

    P=Pext*760

    x = lambda T: (P-Psat[FluidoB](T))/(Psat[FluidoA](T)-Psat[FluidoB](T))
    y = lambda T: x(T)*Psat[FluidoA](T)/P
    T = np.linspace(Tsat[FluidoA](P),Tsat[FluidoB](P),100)
    plt.figure(figsize=(5,5),dpi=100)
    xT=x(T)
    yT=y(T)
    p=xT.argsort()
    xT=xT[p]
    yT=yT[p]
    plt.plot(xT,yT, color='black')
    plt.plot([0,1],[0,1],color='black',linestyle='--')
    plt.axis('equal')

    # RECTA Q
    Tdew = fsolve(lambda T: y(T)-xF, 138)[0]
    Tbub = fsolve(lambda T: x(T)-xF, 0.01)[0]
    
    Hvapmix=xF*Hvap[FluidoA]+(1-xF)*Hvap[FluidoB]
    cPVapmix=xF*cPVap[FluidoA]+(1-xF)*cPVap[FluidoB]
    cPLiqmix=xF*cPLiq[FluidoA]+(1-xF)*cPLiq[FluidoB]
    
    if (TF>Tdew):
        q = -(TF-Tdew)*(cPVapmix)/Hvapmix
    elif (TF>Tbub):
        q = (Tdew-TF)/(Tdew-Tbub)
    else:
        q = (Hvapmix+(Tbub-TF)*cPLiqmix)/Hvapmix
    
    yq=-xT*q/(1-q)+xF/(1-q)

    #Intersection occurs between points itemindex and itemindex+1
    itemindex = np.argwhere(np.diff(np.sign(yq - yT))).flatten()
    a1=np.array([xT[itemindex][0],yq[itemindex][0]])
    a2=np.array([xT[itemindex+1][0],yq[itemindex+1][0]])
    b1=np.array([xT[itemindex][0],yT[itemindex][0]])
    b2=np.array([xT[itemindex+1][0],yT[itemindex+1][0]])
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = np.empty_like(da)
    dap[0] = -da[1]
    dap[1] = da[0]
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    intersec=(num / denom.astype(float))*db + b1    
    plt.plot([xF,intersec[0]],[xF,intersec[1]],color='blue')
    plt.plot(intersec[0],intersec[1],'ko',ms=5)

    # RECTA DE OPERACION DE ALIMENTACION (REFLUJO MINIMO)
    aux= xD-xD*(xD-intersec[1])/(xD-intersec[0])
    
    #plt.plot([xD, intersec[0]], [xD,intersec[1]], 'c--')
    plt.plot([xD, 0], [xD,aux], 'c--')
    Rextmin=(intersec[1]-xD)/(intersec[0]-intersec[1])
    
    if (Rext<Rextmin):
        print('Reflux less than the minimum allowed (%.4f)'%Rextmin)
    
    # RECTA DE OPERACION DE ALIMENTACION
    plt.plot([xD,0],[xD,xD/(1+Rext)],'c-')
    
    # INTERSECCION RECTA DE OPERACION - RECTA DE ALIMENTACION
    a1=np.array([xF,xF])
    a2=intersec
    b1=np.array([xD,xD])
    b2=np.array([0,xD/(1+Rext)])
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = np.empty_like(da)
    dap[0] = -da[1]
    dap[1] = da[0]
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    intersec2=(num / denom.astype(float))*db + b1    
    plt.plot(intersec2[0],intersec2[1],'mo',ms=5)
    
    plt.title('x-y diagram for {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
    plt.xlabel("$x_{"+FluidoA+"}$")
    plt.ylabel("$y_{"+FluidoA+"}$")

    plt.plot(xF,xF,'go',ms=5)
    plt.plot(xD,xD,'bo',ms=5)

    plt.annotate('$x_D/(1+R_{ext})$', xy=(0.05, xD/(1+Rext)/2), xytext=(0.1, xD/(1+Rext)/2), 
            ha='left', va='center',
            arrowprops=dict(arrowstyle='-[, widthB=%f, lengthB=1.5'%(14.0*xD/(1+Rext)), lw=1.0, color='g'), color='g')    
    
    plt.annotate('$x_D/(1+R_{ext,min})$\n$R_{ext,min}=%.2f$'%Rextmin, xy=(0.35, aux/2), xytext=(0.45, aux/2), 
            ha='left', va='center',
            arrowprops=dict(arrowstyle='-[, widthB=%f, lengthB=10.0'%(14.0*aux), lw=1.0, color='r'), color='r')    

    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.minorticks_on()
    plt.grid(linewidth=1, which='both')

    
widgetFeedRectif=interactive(checkFeedRectif,
                         TF=widgets.BoundedFloatText(value=100,min=30,max=300,step=5,description='$T_F$ (°C)'), 
                         xF=widgets.BoundedFloatText(value=0.5,min=0.1,max=0.9,step=0.01,description='$x_F$'), 
                         xD=widgets.BoundedFloatText(value=0.9,min=0.51,max=0.99,step=0.01,description='$x_D$'),                           Rext=widgets.BoundedFloatText(value=3.0,min=0.1,max=5.0,step=0.1,description='$R_{ext}$'))
controlsFeedRectif = HBox(widgetFeedRectif.children[:-1], layout = Layout(flex_flow='row wrap'))
outputFeedRectif = widgetFeedRectif.children[-1]

##### Diagrama y-x con todas las líneas de operación necearias

def checkAllLines(TF=100.0,xF=0.5,xB=0.1,xD=0.9,Rext=3.0):
    FluidoA='Isopropyl-alcohol'
    FluidoB='Acetic acid'
    Pext=2.0

    P=Pext*760
    x = lambda T: (P-Psat[FluidoB](T))/(Psat[FluidoA](T)-Psat[FluidoB](T))
    y = lambda T: x(T)*Psat[FluidoA](T)/P
    T = np.linspace(Tsat[FluidoA](P),Tsat[FluidoB](P),100)
    plt.figure(figsize=(5,5),dpi=100)
    xT=x(T)
    yT=y(T)
    p=xT.argsort()
    xT=xT[p]
    yT=yT[p]
    plt.plot(xT,yT, color='black')
    plt.plot([0,1],[0,1],color='black',linestyle='--')
    plt.axis('equal')

    # RECTA Q
    Tdew = fsolve(lambda T: y(T)-xF, 138)[0]
    Tbub = fsolve(lambda T: x(T)-xF, 0.01)[0]
    
    Hvapmix=xF*Hvap[FluidoA]+(1-xF)*Hvap[FluidoB]
    cPVapmix=xF*cPVap[FluidoA]+(1-xF)*cPVap[FluidoB]
    cPLiqmix=xF*cPLiq[FluidoA]+(1-xF)*cPLiq[FluidoB]
    
    if (TF>Tdew):
        q = -(TF-Tdew)*(cPVapmix)/Hvapmix
    elif (TF>Tbub):
        q = (Tdew-TF)/(Tdew-Tbub)
    else:
        q = (Hvapmix+(Tbub-TF)*cPLiqmix)/Hvapmix
    
    yq=-xT*q/(1-q)+xF/(1-q)

    #Intersection occurs between points itemindex and itemindex+1
    itemindex = np.argwhere(np.diff(np.sign(yq - yT))).flatten()
    a1=np.array([xT[itemindex][0],yq[itemindex][0]])
    a2=np.array([xT[itemindex+1][0],yq[itemindex+1][0]])
    b1=np.array([xT[itemindex][0],yT[itemindex][0]])
    b2=np.array([xT[itemindex+1][0],yT[itemindex+1][0]])
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = np.empty_like(da)
    dap[0] = -da[1]
    dap[1] = da[0]
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    intersec=(num / denom.astype(float))*db + b1    
    #plt.plot([xF,intersec[0]],[xF,intersec[1]],color='blue')
    #plt.plot(intersec[0],intersec[1],'ko',ms=5)

    # RECTA DE OPERACION DE ALIMENTACION (REFLUJO MINIMO)
    #plt.plot([xD, intersec[0]], [xD,intersec[1]], 'c--')
    Rextmin=(intersec[1]-xD)/(intersec[0]-intersec[1])
    
    permitted=True
    if (Rext<Rextmin):
        print('Reflux less than the minimum allowed (%.4f)'%Rextmin)
        permitted=False
        
    # INTERSECCION RECTA DE OPERACION - RECTA DE ALIMENTACION
    a1=np.array([xF,xF])
    a2=intersec
    b1=np.array([xD,xD])
    b2=np.array([0,xD/(1+Rext)])
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = np.empty_like(da)
    dap[0] = -da[1]
    dap[1] = da[0]
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    intersec2=(num / denom.astype(float))*db + b1    
    if permitted:
        plt.plot(intersec2[0],intersec2[1],'mo',ms=5)
        plt.plot([xF,intersec2[0]],[xF,intersec2[1]],color='blue')
        # RECTA DE OPERACION DE ALIMENTACION
        plt.plot([xD,intersec2[0]],[xD,intersec2[1]],'c-')
    
        # RECTA DE OPERACION - ZONA DE AGOTAMIENTO
        plt.plot([xB,intersec2[0]],[xB,intersec2[1]])
    else:
        plt.plot(intersec2[0],intersec2[1],'ro',ms=5)
        plt.plot([xF,intersec2[0]],[xF,intersec2[1]],color='r')
        # RECTA DE OPERACION DE ALIMENTACION
        plt.plot([xD,intersec2[0]],[xD,intersec2[1]],'r-')
    
        # RECTA DE OPERACION - ZONA DE AGOTAMIENTO
        plt.plot([xB,intersec2[0]],[xB,intersec2[1]],'r-')
        
        
    plt.title('x-y diagram for {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
    plt.xlabel("$x_{"+FluidoA+"}$")
    plt.ylabel("$y_{"+FluidoA+"}$")

    plt.plot(xF,xF,'go',ms=5)
    plt.plot(xD,xD,'bo',ms=5)
    plt.plot(xB,xB,'yo',ms=5)

    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.minorticks_on()
    plt.grid(linewidth=1, which='both')

    
widgetAllLines=interactive(checkAllLines,
                         TF=widgets.BoundedFloatText(value=100,min=30,max=300,step=5,description='$T_F$ (°C)'), 
                         xF=widgets.BoundedFloatText(value=0.5,min=0.1,max=0.9,step=0.01,description='$x_F$'), 
                         xB=widgets.BoundedFloatText(value=0.1,min=0.01,max=0.49,step=0.01,description='$x_B$'), 
                         xD=widgets.BoundedFloatText(value=0.9,min=0.51,max=0.99,step=0.01,description='$x_D$'),                           Rext=widgets.BoundedFloatText(value=3.0,min=0.1,max=5.0,step=0.1,description='$R_{ext}$'))
controlsAllLines = HBox(widgetAllLines.children[:-1], layout = Layout(flex_flow='row wrap'))
outputAllLines = widgetAllLines.children[-1]

##### Diagrama completo McCabe-Thiele paso a paso

def checkMcCabeStep(TF=100.0,xF=0.5,xB=0.1,xD=0.9,Rext=3.0,Steps=1):

    FluidoA='Isopropyl-alcohol'
    FluidoB='Acetic acid'
    Pext=2.0

    P=Pext*760
    x = lambda T: (P-Psat[FluidoB](T))/(Psat[FluidoA](T)-Psat[FluidoB](T))
    y = lambda T: x(T)*Psat[FluidoA](T)/P
    T = np.linspace(Tsat[FluidoA](P),Tsat[FluidoB](P),100)
    #plt.figure(figsize=(5,5),dpi=100)
    xT=x(T)
    yT=y(T)
    p=xT.argsort()
    xT=xT[p]
    yT=yT[p]
    T=T[p]

    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_axes((.05,.1,.45,.8))    # coord x, coord y, ancho, alto
    plt.plot(xT,yT, color='black')
    plt.plot([0,1],[0,1],color='black',linestyle='--')
    plt.axis('equal')

    # RECTA Q
    Tdew = fsolve(lambda T: y(T)-xF, 138)[0]
    Tbub = fsolve(lambda T: x(T)-xF, 0.01)[0]

    Hvapmix=xF*Hvap[FluidoA]+(1-xF)*Hvap[FluidoB]
    cPVapmix=xF*cPVap[FluidoA]+(1-xF)*cPVap[FluidoB]
    cPLiqmix=xF*cPLiq[FluidoA]+(1-xF)*cPLiq[FluidoB]

    if (TF>Tdew):
        q = -(TF-Tdew)*(cPVapmix)/Hvapmix
    elif (TF>Tbub):
        q = (Tdew-TF)/(Tdew-Tbub)
    else:
        q = (Hvapmix+(Tbub-TF)*cPLiqmix)/Hvapmix

    yq=-xT*q/(1-q)+xF/(1-q)

    #Intersection occurs between points itemindex and itemindex+1
    itemindex = np.argwhere(np.diff(np.sign(yq - yT))).flatten()
    a1=np.array([xT[itemindex][0],yq[itemindex][0]])
    a2=np.array([xT[itemindex+1][0],yq[itemindex+1][0]])
    b1=np.array([xT[itemindex][0],yT[itemindex][0]])
    b2=np.array([xT[itemindex+1][0],yT[itemindex+1][0]])
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = np.empty_like(da)
    dap[0] = -da[1]
    dap[1] = da[0]
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    intersec=(num / denom.astype(float))*db + b1    
    #plt.plot([xF,intersec[0]],[xF,intersec[1]],color='blue')
    #plt.plot(intersec[0],intersec[1],'ko',ms=5)

    # RECTA DE OPERACION DE ALIMENTACION (REFLUJO MINIMO)
    #plt.plot([xD, intersec[0]], [xD,intersec[1]], 'c--')
    Rextmin=(intersec[1]-xD)/(intersec[0]-intersec[1])

    permitted=True
    if (Rext<Rextmin):
        print('Reflux less than the minimum allowed (%.4f)'%Rextmin)
        permitted=False

    # INTERSECCION RECTA DE OPERACION - RECTA DE ALIMENTACION
    a1=np.array([xF,xF])
    a2=intersec
    b1=np.array([xD,xD])
    b2=np.array([0,xD/(1+Rext)])
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = np.empty_like(da)
    dap[0] = -da[1]
    dap[1] = da[0]
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    intersec2=(num / denom.astype(float))*db + b1    

    plt.title('x-y diagram for \n{:s}/{:s} at p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
    plt.xlabel("$x_{"+FluidoA+"}$")
    plt.ylabel("$y_{"+FluidoA+"}$")

    plt.plot(xF,xF,'go',ms=5)
    plt.plot(xD,xD,'bo',ms=5)
    plt.plot(xB,xB,'yo',ms=5)

    if permitted:
        plt.plot(intersec2[0],intersec2[1],'mo',ms=5)
        plt.plot([xF,intersec2[0]],[xF,intersec2[1]],color='blue')
        # RECTA DE OPERACION DE ALIMENTACION
        plt.plot([xD,intersec2[0]],[xD,intersec2[1]],'c-')

        # RECTA DE OPERACION - ZONA DE AGOTAMIENTO
        plt.plot([xB,intersec2[0]],[xB,intersec2[1]])
    else:
        plt.plot(intersec2[0],intersec2[1],'ro',ms=5)
        plt.plot([xF,intersec2[0]],[xF,intersec2[1]],color='r')
        # RECTA DE OPERACION DE ALIMENTACION
        plt.plot([xD,intersec2[0]],[xD,intersec2[1]],'r-')

        # RECTA DE OPERACION - ZONA DE AGOTAMIENTO
        plt.plot([xB,intersec2[0]],[xB,intersec2[1]],'r-')


    S = (intersec2[1]-xB)/(intersec2[0]-xB)

    # START McCabe-Thiele
    xP = xD
    yP = xD
    nTray = 0

    while xP > xB and nTray < Steps:
        nTray += 1
        Tdew = fsolve(lambda T:y(T) - yP, 100)
        xQ = xP
        xP = float(x(Tdew))
        plt.plot([xQ,xP],[yP,yP],'r')
        #plt.plot(xP,yP,'ro',ms=5)
        plt.text(xP-0.03,yP,nTray)

        yQ = yP
        yP = min([xD - (Rext/(Rext+1))*(xD-xP),xB + S*(xP-xB)])
        plt.plot([xP,xP],[yQ,yP],'r')

    plt.plot([xB,xB],[0,xB],'b--')
    plt.plot([xF,xF],[0,xF],'b--')
    plt.plot([xD,xD],[0,xD],'b--')
    plt.text(xB,0.02,'$x_B$ = {:0.2f}'.format(float(xB)),horizontalalignment='center')
    plt.text(xF,0.02,'$x_F$ = {:0.2f}'.format(float(xF)),horizontalalignment='center')
    plt.text(xD,0.02,'$x_D$ = {:0.2f}'.format(float(xD)),horizontalalignment='center')

    if xP <= xB:
        plt.text(0.05,0.9,'{:d} steps eq.\nRmin={:.4g}'.format(int(nTray),Rextmin),weight='bold')
    else:
        plt.text(0.05,0.9,'\nRmin={:.4g}'.format(Rextmin),weight='bold')

    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.minorticks_on()
    plt.grid(linewidth=1, which='both')

    ########################################################################################################################
    #  fin codigo original, adicion de nuevo grafico. He modificado el formato de grafico a figure y distribucion con axes #
    ########################################################################################################################

    # xT, yT y T son nparrays de 100 puntos. Necesito las posiciones de los platos en el vector de 100 puntos

    #return xT, yT

    # repeticion del bucle de los platos, en el futuro se puede borrar e introducir en el original

    # 2 POSOBILIDADES: Interpolacion o representacion directa

    caso="interpol"
    caso="Repre directa"

    x_plato=[]
    y_plato=[]
    T_plato=[]
    isFeed = -1

    if caso=="interpol":
        xP = xD
        yP = xD
        nTray = 0

        x_plato.append(xP)    # punto de equilibrio en plato 0

        while xP > xB:
            nTray += 1
            Tdew = fsolve(lambda T:y(T) - yP, 100)
            xQ = xP
            xP = x(Tdew)
            #plt.plot([xQ,xP],[yP,yP],'r')
            #plt.plot(xP,yP,'ro',ms=5)
            plt.text(xP-0.03,yP,nTray)

            yQ = yP
            yP1 = xD - (Rext/(Rext+1))*(xD-xP)
            yP2 = xB + S*(xP-xB)
            if yP1 < yP2:
                yP = yP1
            else:
                yP = yP2
                if (isFeed<0):
                    isFeed=nTray
            #yP = min([xD - (Rext/(Rext+1))*(xD-xP),xB + S*(xP-xB)])
            #plt.plot([xP,xP],[yQ,yP],'r')

            x_plato.append(xP)
        x_plato.append(x(fsolve(lambda T:y(T) - yP, 100)))   # punto de equilibirio en plato posterior al final


        comprobacion=[]          # para ver si el valor de posicion es igual al de los platos
        posicion=[]

            # bucle para posiciones de plato dentro del vector

        for i in x_plato:        # bucle a lo largo de los platos
                n=0
                for j in range(0,len(xT)):        # bucle a lo largo de xT 
                    if xT[j]>i and n==0:
                        comprobacion.append(xT[j])
                        posicion.append(j)    # esta es la posicion donde xT es mayor que la x calculada
                        n=1
        for i in posicion:
                y_plato.append(yT[i])
                T_plato.append(T[i])



    if caso=="Repre directa":
        xP = xD
        yP = xD
        nTray = 0

        #T_plato_0=fsolve(lambda T: xD*Psat[FluidoA](T)/P, 100)# aqui no se calcular la temperatura del plato 0. T cuando x=xD
        T_plato_0=fsolve(lambda T: x(T)-xP,100)
        x_plato.append(xD)          # x en el plato 0
        y0=y(T_plato_0)
        y_plato.append(float(y0))        # y en el plato 0
        T_plato.append(float(T_plato_0))   #T en el plato 0

            # la composicion del liquido en cada plato es x(Tdew)
        extra=2
        while (xP > xB) or extra > 0:  # calculo para todos los platos + 1
            nTray += 1
            Tdew = fsolve(lambda T:y(T) - yP, 100)
            T_plato.append(float(Tdew))
            xQ = xP
            xP = x(Tdew)

            #plt.plot([xQ,xP],[yP,yP],'r')
            #plt.plot(xP,yP,'ro',ms=5)
            #plt.text(xP-0.03,yP,nTray)

            yQ = yP
            yP1 = xD - (Rext/(Rext+1))*(xD-xP)
            yP2 = xB + S*(xP-xB)
            if yP1 < yP2:
                yP = yP1
            else:
                yP = yP2
                if (isFeed<0):
                    isFeed=nTray
            # yP = min([xD - (Rext/(Rext+1))*(xD-xP),xB + S*(xP-xB)])
            #plt.plot([xP,xP],[yQ,yP],'r')

            x_plato.append(float(xP))
            y_plato.append(float(yQ))

            if xP<=xB:
                extra -=1   #composicion del liquido y vapor si se hiciese un plato mas

        nTray -=1        

        #plt.plot([xB,xB],[0,xB],'b--')
        #plt.plot([xF,xF],[0,xF],'b--')
        #plt.plot([xD,xD],[0,xD],'b--')
        #plt.text(xB,0.02,'$x_B$ = {:0.2f}'.format(float(xB)),horizontalalignment='center')
        #plt.text(xF,0.02,'$x_F$ = {:0.2f}'.format(float(xF)),horizontalalignment='center')
        #plt.text(xD,0.02,'$x_D$ = {:0.2f}'.format(float(xD)),horizontalalignment='center')








    # valores para el zoom
    if Steps<=nTray:
        n0=Steps-1
        n1=Steps
        n2=Steps+1
    else:
        n0=nTray-1
        n1=nTray
        n2=nTray+1

    if x_plato[n1]<0.1:
        eje=[0, y_plato[n1]+0.1,T_plato[n0]-(T_plato[n2]-T_plato[n1])/2, T_plato[n2]+(T_plato[n1]-T_plato[n0])/2]
    elif x_plato[n0]>0.9:
        eje=[x_plato[n1]-0.1, 1,T_plato[n2]-(T_plato[n0]-T_plato[n1])/2, T_plato[n0]+(T_plato[n0]-T_plato[n1])/2]
    else:
        eje=[x_plato[n1]-0.1, y_plato[n1]+0.1,T_plato[n0]-(T_plato[n2]-T_plato[n1])/2, T_plato[n2]+(T_plato[n1]-T_plato[n0])/2]

    #eje=[0,0.5,90,100]



    #representacion del equilibrio (zoom)
    if permitted:
        ax2 = fig.add_axes((.55,.1,.3,.6))
        plt.xlabel('x,y (%s)'%FluidoA)
        plt.ylabel('T')
        plt.title('T-xy diagram\nat plate {:d}'.format(n1),loc='left')
        #plt.plot(xT,T,yT,T)
        T=T[p]
        plt.plot([x(T) for T in T],T)
        plt.plot([y(T) for T in T],T)

        plt.plot([x_plato[n1],y_plato[n1]],[T_plato[n1],T_plato[n1]], color='navy', alpha=.75, lw=2, ls='dotted')
        plt.plot([x_plato[n1-1],y_plato[n1+1]],[T_plato[n1-1],T_plato[n1+1]], color='forestgreen',alpha=.75, lw=2, ls='dotted' )
        plt.scatter([x_plato[n1],y_plato[n1],x_plato[n1-1],y_plato[n1+1]],[T_plato[n1],T_plato[n1],T_plato[n1-1],T_plato[n1+1]])
        plt.axis(eje)

            # flechas liquido vapor
        dist_hor=float(abs(eje[0]-eje[1]))
        dist_vert=float(abs(eje[2]-eje[3]))
        dx=dist_hor/4
        dy=dist_vert/4
        ax2.arrow(x_plato[n1]+dist_hor*0.05,T_plato[n1],0,-1*dy,head_width=0.01,head_length=1,color='navy')  #x,y,dx, dy, 
        ax2.arrow(y_plato[n1]-dist_hor*0.05,T_plato[n1],0,dy,head_width=0.01,head_length=1,color='#cc5801')

        ax2.text(x_plato[n1]+dist_hor*0.08,T_plato[n1]-dist_vert*0.25,'L',color='navy')
        ax2.text(y_plato[n1]-dist_hor*0.1,T_plato[n1]+dist_vert*0.25,'V',color='#cc5801')
            # composiciones
        ax2.text(x_plato[n1]-dist_hor*0.05,T_plato[n1]+dist_vert*0.07,("$x_{%d}$ = %.2f"%(n1,round(x_plato[n1],2))))
        ax2.text(y_plato[n1]-dist_hor*0.05,T_plato[n1]-dist_vert*0.1,("$y_{%d}$ = %.2f"%(n1,round(y_plato[n1],2))))
        ax2.text(x_plato[n1-1]+dist_hor*0.05,T_plato[n1-1],("$x_{%d}$ = %.2f"%(n0,round(x_plato[n1-1],2))))
        ax2.text(y_plato[n1+1]+dist_hor*0.05,T_plato[n1+1],("$y_{%d}$ = %.2f"%(n2,round(x_plato[n1+1],2))))
        #axes.legend(loc=0)


        #plt.plot([x_platos[n1],y_platos[n1]],[T_platos[n1],T_platos[n1]], color='navy', alpha=.75, lw=2, ls='dotted')
        #plt.plot([x_platos[n0],y_platos[n2]],[T_platos[n0],T_platos[n2]], color='forestgreen',alpha=.75, lw=2, ls='dotted' )
        #plt.scatter([x_platos[n1],y_platos[n1],x_platos[n0],y_platos[n2]],[T_platos[n1],T_platos[n1],T_platos[n0],T_platos[n2]])
        #plt.axis(eje)

        ax3 = fig.add_axes((0.8,0.6,0.18,0.35))

        plt.plot(xT[p],T,yT[p],T)
        plt.plot([x_plato[n1],y_plato[n1]],[T_plato[n1],T_plato[n1]], color='navy', alpha=.75, lw=2, ls='dotted')
        plt.plot([x_plato[n0],y_plato[n2]],[T_plato[n0],T_plato[n2]], color='forestgreen',alpha=.75, lw=2, ls='dotted' )
        plt.plot([xF,xF],[T[0],TF],color='b',ls='dashed', alpha=0.5)
        plt.plot([0,xF],[TF,TF],color='b',ls='dashed', alpha=0.5)
        plt.scatter(xF,TF, marker='X',alpha=1, edgecolors='yellow')

        info=('$T_{%d}$=%.1f'%(n1,T_plato[n1]))
        plt.plot([], [], ' ', label=info)
        plt.legend()
        ax3.text(xF+0.05,TF,'F',color='navy')
    
        # Mini-column
        # Assume F=1: 
        # Mass balance
        # F = D + B
        # xF*F = xD*D + xB*B = xD*(F-B)+ xB*B = xD*F+B*(xB-xD) --> B = F*(xD-xF)/(xD-xB) --> D = F*(xF-xB)/(xD-xB)
        B = (xD-xF)/(xD-xB)
        D = (xF-xB)/(xD-xB)
        L = Rext*D

        ax4 = fig.add_axes((0.85,0.1,0.15,0.4))
        ax4.axis('off')
        R=L
        H=nTray
        theta=np.linspace(np.pi, 0, 100)
        xbase = R*np.cos(theta)
        ybase = R*np.sin(theta)
        xcyl_right = np.array([R, R])
        ycyl_right = np.array([0, -H])
        theta=np.linspace(0, -np.pi, 100)
        xhead = R*np.cos(theta)
        yhead = -H + R*np.sin(theta)
        xcyl_left = np.array([-R, -R])
        ycyl_left = np.array([-H, 0])
        ax4.plot(np.concatenate((xbase, xcyl_right, xhead, xcyl_left)), np.concatenate((ybase, ycyl_right, yhead, ycyl_left)))
        plt.xlim(-5, 5)
        for i in range(nTray):
            if i == n1-1:
                col='r'
                wid=2
            else:
                col='k'
                wid=1
            ax4.plot([-R*1.1, R*1.1], [-0.5-i, -0.5 - i], color=col, linewidth=wid)
            ax4.text(R*1.2, -0.55-i, "%d"%(i+1), color=col)
        ax4.arrow(-R-0.5,-isFeed+0.75,0.45,0,head_width=0.25,head_length=0.25, color='navy')  #x,y,dx, dy, 
        ax4.text(-R-0.5,-isFeed+0.8,"F")
    
    plt.show()


    
widgetMcCabeStep=interactive(checkMcCabeStep, manual=True,
                         TF=widgets.BoundedFloatText(value=100,min=30,max=300,step=5,description='$T_F$ (°C)'), 
                         xF=widgets.BoundedFloatText(value=0.5,min=0.1,max=0.9,step=0.01,description='$x_F$'), 
                         xB=widgets.BoundedFloatText(value=0.1,min=0.01,max=0.49,step=0.01,description='$x_B$'), 
                         xD=widgets.BoundedFloatText(value=0.9,min=0.51,max=0.99,step=0.01,description='$x_D$'),                           
                         Rext=widgets.BoundedFloatText(value=3.0,min=0.1,max=5.0,step=0.1,description='$R_{ext}$'),
                         Steps=widgets.BoundedIntText(value=1,min=1,max=50,step=1))
controlsMcCabeStep = HBox(widgetMcCabeStep.children[:-1], layout = Layout(flex_flow='row wrap'))
outputMcCabeStep = widgetMcCabeStep.children[-1]

##### Diagrama completo McCabe-Thiele

def checkMcCabe(FluidoA='Acetone', FluidoB='Water',Pext=1.0,TF=100.0,xF=0.5,xB=0.1,xD=0.9,Rext=3.0):
    if FluidoA==FluidoB:
        print("Choose two different fluids")
        return
    P=Pext*760
    x = lambda T: (P-Psat[FluidoB](T))/(Psat[FluidoA](T)-Psat[FluidoB](T))
    y = lambda T: x(T)*Psat[FluidoA](T)/P
    T = np.linspace(Tsat[FluidoA](P),Tsat[FluidoB](P),100)
    plt.figure(figsize=(5,5),dpi=100)
    xT=x(T)
    yT=y(T)
    p=xT.argsort()
    xT=xT[p]
    yT=yT[p]
    plt.plot(xT,yT, color='black')
    plt.plot([0,1],[0,1],color='black',linestyle='--')
    plt.axis('equal')
    if yT[50]<xT[50]:
        print("Select FluidA as the most volatile fluid")
        return

    # RECTA Q
    Tdew = fsolve(lambda T: y(T)-xF, 138)[0]
    Tbub = fsolve(lambda T: x(T)-xF, 0.01)[0]
    
    Hvapmix=xF*Hvap[FluidoA]+(1-xF)*Hvap[FluidoB]
    cPVapmix=xF*cPVap[FluidoA]+(1-xF)*cPVap[FluidoB]
    cPLiqmix=xF*cPLiq[FluidoA]+(1-xF)*cPLiq[FluidoB]
    
    if (TF>Tdew):
        q = -(TF-Tdew)*(cPVapmix)/Hvapmix
    elif (TF>Tbub):
        q = (Tdew-TF)/(Tdew-Tbub)
    else:
        q = (Hvapmix+(Tbub-TF)*cPLiqmix)/Hvapmix
    
    yq=-xT*q/(1-q)+xF/(1-q)

    #Intersection occurs between points itemindex and itemindex+1
    itemindex = np.argwhere(np.diff(np.sign(yq - yT))).flatten()
    a1=np.array([xT[itemindex][0],yq[itemindex][0]])
    a2=np.array([xT[itemindex+1][0],yq[itemindex+1][0]])
    b1=np.array([xT[itemindex][0],yT[itemindex][0]])
    b2=np.array([xT[itemindex+1][0],yT[itemindex+1][0]])
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = np.empty_like(da)
    dap[0] = -da[1]
    dap[1] = da[0]
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    intersec=(num / denom.astype(float))*db + b1    
    #plt.plot([xF,intersec[0]],[xF,intersec[1]],color='blue')
    #plt.plot(intersec[0],intersec[1],'ko',ms=5)

    # RECTA DE OPERACION DE ALIMENTACION (REFLUJO MINIMO)
    #plt.plot([xD, intersec[0]], [xD,intersec[1]], 'c--')
    Rextmin=(intersec[1]-xD)/(intersec[0]-intersec[1])
    
    permitted=True
    if (Rext<Rextmin):
        print('Reflux less than the minimum allowed (%.4f)'%Rextmin)
        permitted=False
        
    # INTERSECCION RECTA DE OPERACION - RECTA DE ALIMENTACION
    a1=np.array([xF,xF])
    a2=intersec
    b1=np.array([xD,xD])
    b2=np.array([0,xD/(1+Rext)])
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = np.empty_like(da)
    dap[0] = -da[1]
    dap[1] = da[0]
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    intersec2=(num / denom.astype(float))*db + b1    
    
    plt.title('x-y diagram for {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
    plt.xlabel("$x_{"+FluidoA+"}$")
    plt.ylabel("$y_{"+FluidoA+"}$")

    plt.plot(xF,xF,'go',ms=5)
    plt.plot(xD,xD,'bo',ms=5)
    plt.plot(xB,xB,'yo',ms=5)

    if permitted:
        plt.plot(intersec2[0],intersec2[1],'mo',ms=5)
        plt.plot([xF,intersec2[0]],[xF,intersec2[1]],color='blue')
        # RECTA DE OPERACION DE ALIMENTACION
        plt.plot([xD,intersec2[0]],[xD,intersec2[1]],'c-')
    
        # RECTA DE OPERACION - ZONA DE AGOTAMIENTO
        plt.plot([xB,intersec2[0]],[xB,intersec2[1]])
    else:
        plt.plot(intersec2[0],intersec2[1],'ro',ms=5)
        plt.plot([xF,intersec2[0]],[xF,intersec2[1]],color='r')
        # RECTA DE OPERACION DE ALIMENTACION
        plt.plot([xD,intersec2[0]],[xD,intersec2[1]],'r-')
    
        # RECTA DE OPERACION - ZONA DE AGOTAMIENTO
        plt.plot([xB,intersec2[0]],[xB,intersec2[1]],'r-')
        return        
        
    S = (intersec2[1]-xB)/(intersec2[0]-xB)

    # START McCabe-Thiele
    xP = xD
    yP = xD
    nTray = 0

    while xP > xB:
        nTray += 1
        Tdew = fsolve(lambda T:y(T) - yP, 100)
        xQ = xP
        xP = x(Tdew)
        plt.plot([xQ,xP],[yP,yP],'r')
        #plt.plot(xP,yP,'ro',ms=5)
        plt.text(xP-0.03,yP,nTray)

        yQ = yP
        yP = min([xD - (Rext/(Rext+1))*(xD-xP),xB + S*(xP-xB)])
        plt.plot([xP,xP],[yQ,yP],'r')
        
    plt.plot([xB,xB],[0,xB],'b--')
    plt.plot([xF,xF],[0,xF],'b--')
    plt.plot([xD,xD],[0,xD],'b--')
    plt.text(xB,0.02,'$x_B$ = {:0.2f}'.format(float(xB)),horizontalalignment='center')
    plt.text(xF,0.02,'$x_F$ = {:0.2f}'.format(float(xF)),horizontalalignment='center')
    plt.text(xD,0.02,'$x_D$ = {:0.2f}'.format(float(xD)),horizontalalignment='center')
    
    plt.text(0.05,0.9,'{:d} steps eq.\nRmin={:.4g}'.format(int(nTray),Rextmin),weight='bold')

    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.minorticks_on()
    plt.grid(linewidth=1, which='both')

    
widgetMcCabe=interactive(checkMcCabe,FluidoA=Psat.keys(),FluidoB=Psat.keys(),
                         Pext=widgets.BoundedFloatText(value=1.0,min=0.1,max=5.0,step=0.1,description='p (atm)'), 
                         TF=widgets.BoundedFloatText(value=100,min=30,max=300,step=5,description='$T_F$ (°C)'), 
                         xF=widgets.BoundedFloatText(value=0.5,min=0.1,max=0.9,step=0.01,description='$x_F$'), 
                         xB=widgets.BoundedFloatText(value=0.1,min=0.01,max=0.49,step=0.01,description='$x_B$'), 
                         xD=widgets.BoundedFloatText(value=0.9,min=0.51,max=0.99,step=0.01,description='$x_D$'),                           Rext=widgets.BoundedFloatText(value=3.0,min=0.1,max=5.0,step=0.1,description='$R_{ext}$'))
controlsMcCabe = HBox(widgetMcCabe.children[:-1], layout = Layout(flex_flow='row wrap'))
outputMcCabe = widgetMcCabe.children[-1]

##### Diagrama incompleto McCabe-Thiele

def checkMcCabeIncomplete(FluidoA='Acetone', FluidoB='Water',Pext=1.0,TF=100.0,xF=0.5,xB=0.1,xD=0.9,Rext=3.0):
    if FluidoA==FluidoB:
        print("Choose two different fluids")
        return
    P=Pext*760
    x = lambda T: (P-Psat[FluidoB](T))/(Psat[FluidoA](T)-Psat[FluidoB](T))
    y = lambda T: x(T)*Psat[FluidoA](T)/P
    T = np.linspace(Tsat[FluidoA](P),Tsat[FluidoB](P),100)
    plt.figure(figsize=(5,5),dpi=300)
    xT=x(T)
    yT=y(T)
    p=xT.argsort()
    xT=xT[p]
    yT=yT[p]
    plt.plot(xT,yT, color='black')
    plt.plot([0,1],[0,1],color='black',linestyle='--')
    plt.axis('equal')
    if yT[50]<xT[50]:
        print("Select FluidA as the most volatile fluid")
        return

    # RECTA Q
    Tdew = fsolve(lambda T: y(T)-xF, 138)[0]
    Tbub = fsolve(lambda T: x(T)-xF, 0.01)[0]
    
    Hvapmix=xF*Hvap[FluidoA]+(1-xF)*Hvap[FluidoB]
    cPVapmix=xF*cPVap[FluidoA]+(1-xF)*cPVap[FluidoB]
    cPLiqmix=xF*cPLiq[FluidoA]+(1-xF)*cPLiq[FluidoB]
    
    if (TF>Tdew):
        q = -(TF-Tdew)*(cPVapmix)/Hvapmix
    elif (TF>Tbub):
        q = (Tdew-TF)/(Tdew-Tbub)
    else:
        q = (Hvapmix+(Tbub-TF)*cPLiqmix)/Hvapmix
    
    yq=-xT*q/(1-q)+xF/(1-q)

    #Intersection occurs between points itemindex and itemindex+1
    itemindex = np.argwhere(np.diff(np.sign(yq - yT))).flatten()
    a1=np.array([xT[itemindex][0],yq[itemindex][0]])
    a2=np.array([xT[itemindex+1][0],yq[itemindex+1][0]])
    b1=np.array([xT[itemindex][0],yT[itemindex][0]])
    b2=np.array([xT[itemindex+1][0],yT[itemindex+1][0]])
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = np.empty_like(da)
    dap[0] = -da[1]
    dap[1] = da[0]
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    intersec=(num / denom.astype(float))*db + b1    
    plt.plot([xF,intersec[0]],[xF,intersec[1]],color='blue')
    plt.plot(intersec[0],intersec[1],'ko',ms=5)

    # RECTA DE OPERACION DE ALIMENTACION (REFLUJO MINIMO)
    aux= xD-xD*(xD-intersec[1])/(xD-intersec[0])
    
    #plt.plot([xD, intersec[0]], [xD,intersec[1]], 'c--')
    plt.plot([xD, 0], [xD,aux], 'c--')
    Rextmin=(intersec[1]-xD)/(intersec[0]-intersec[1])

    # RECTA DE OPERACION DE ALIMENTACION (REFLUJO MINIMO)
    #plt.plot([xD, intersec[0]], [xD,intersec[1]], 'c--')
    #Rextmin=(intersec[1]-xD)/(intersec[0]-intersec[1])
    
    permitted=True
    if (Rext<Rextmin):
        print('Reflux less than the minimum allowed (%.4f)'%Rextmin)
        permitted=False
    
    # RECTA DE OPERACION DE RECTIFICACION
    plt.plot([xD,0],[xD,xD/(1+Rext)],'c-')

    # INTERSECCION RECTA DE OPERACION - RECTA DE ALIMENTACION
    a1=np.array([xF,xF])
    a2=intersec
    b1=np.array([xD,xD])
    b2=np.array([0,xD/(1+Rext)])
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = np.empty_like(da)
    dap[0] = -da[1]
    dap[1] = da[0]
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    intersec2=(num / denom.astype(float))*db + b1    
    
    plt.title('x-y diagram for {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
    plt.xlabel("$x_{"+FluidoA+"}$")
    plt.ylabel("$y_{"+FluidoA+"}$")

    plt.plot(xF,xF,'go',ms=5)
    plt.plot(xD,xD,'bo',ms=5)
    plt.plot(xB,xB,'yo',ms=5)

    if permitted:
        plt.plot(intersec2[0],intersec2[1],'mo',ms=5)
        plt.plot([xF,intersec2[0]],[xF,intersec2[1]],color='blue')
        # RECTA DE OPERACION DE ALIMENTACION
        plt.plot([xD,intersec2[0]],[xD,intersec2[1]],'c-')
    
        # RECTA DE OPERACION - ZONA DE AGOTAMIENTO
        plt.plot([xB,intersec2[0]],[xB,intersec2[1]])
    else:
        plt.plot(intersec2[0],intersec2[1],'ro',ms=5)
        plt.plot([xF,intersec2[0]],[xF,intersec2[1]],color='r')
        # RECTA DE OPERACION DE ALIMENTACION
        plt.plot([xD,intersec2[0]],[xD,intersec2[1]],'r-')
    
        # RECTA DE OPERACION - ZONA DE AGOTAMIENTO
        plt.plot([xB,intersec2[0]],[xB,intersec2[1]],'r-')
        return        
        
    S = (intersec2[1]-xB)/(intersec2[0]-xB)

    # START McCabe-Thiele
    xP = xD
    yP = xD
    nTray = 0

    while xP > xB:
        nTray += 1
        Tdew = fsolve(lambda T:y(T) - yP, 100)
        xQ = xP
        xP = x(Tdew)
        plt.plot([xQ,xP],[yP,yP],'r')
        #plt.plot(xP,yP,'ro',ms=5)
        #plt.text(xP-0.03,yP,nTray)

        yQ = yP
        yP = min([xD - (Rext/(Rext+1))*(xD-xP),xB + S*(xP-xB)])
        plt.plot([xP,xP],[yQ,yP],'r')
        
    plt.plot([xB,xB],[0,xB],'b--')
    plt.plot([xF,xF],[0,xF],'b--')
    plt.plot([xD,xD],[0,xD],'b--')
    #plt.text(xB,0.02,'$x_B$ = {:0.2f}'.format(float(xB)),horizontalalignment='center')
    #plt.text(xF,0.02,'$x_F$ = {:0.2f}'.format(float(xF)),horizontalalignment='center')
    #plt.text(xD,0.02,'$x_D$ = {:0.2f}'.format(float(xD)),horizontalalignment='center')
    
    #plt.text(0.05,0.9,'{:d} etapas eq.\nRmin={:.4g}'.format(int(nTray),Rextmin),weight='bold')

    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.minorticks_on()
    plt.grid(linewidth=1, which='both')

    
widgetMcCabeIncomplete=interactive(checkMcCabeIncomplete,FluidoA=Psat.keys(),FluidoB=Psat.keys(),
                         Pext=widgets.BoundedFloatText(value=1.0,min=0.1,max=5.0,step=0.1,description='p (atm)'), 
                         TF=widgets.BoundedFloatText(value=100,min=30,max=300,step=5,description='$T_F$ (°C)'), 
                         xF=widgets.BoundedFloatText(value=0.5,min=0.1,max=0.9,step=0.01,description='$x_F$'), 
                         xB=widgets.BoundedFloatText(value=0.1,min=0.01,max=0.49,step=0.01,description='$x_B$'), 
                         xD=widgets.BoundedFloatText(value=0.9,min=0.51,max=0.99,step=0.01,description='$x_D$'),                          Rext=widgets.BoundedFloatText(value=3.0,min=0.1,max=5.0,step=0.1,description='$R_{ext}$'))
controlsMcCabeIncomplete = HBox(widgetMcCabeIncomplete.children[:-1], layout = Layout(flex_flow='row wrap'))
outputMcCabeIncomplete = widgetMcCabeIncomplete.children[-1]
