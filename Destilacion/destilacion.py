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
Psat['Acetona'] = lambda T: 10**(7.1327 - 1219.97/(T + 230.653))
Psat['Acetonitrilo'] = lambda T: 10**(7.33986 - 1482.29/(T + 250.523))
Psat['Acido Acetico'] = lambda T: 10**(7.2996 - 1479.02/(T + 216.82))
Psat['Agua'] = lambda T: 10**(8.07131 - 1730.63/(T + 233.426))
Psat['Etanol'] = lambda T: 10**( 8.20417 - 1642.89/(T + 230.3))
#Psat['Etilenglicol'] = lambda T: 10**( 8.7945 - 2615.4/(T + 244.91))
Psat['Fenol'] = lambda T: 10**( 7.1345 - 1516.07/(T + 174.57))
Psat['isopropyl-alcohol'] = lambda T: 10**( 8.1182 - 1580.92/(T + 219.62))

Molecule = dict()
Molecule['Acetona'] = 'img/acetone.png'
Molecule['Acetonitrilo'] = 'img/acetonitrile.png'
Molecule['Acido Acetico'] = 'img/aceticacid.png'
Molecule['Agua'] = 'img/water.png'
Molecule['Etanol'] = 'img/ethanol.png'
#Molecule['Etilenglicol'] = 'img/ethylenglycol.png'
Molecule['Fenol'] = 'img/phenol.png'
Molecule['isopropyl-alcohol'] = 'img/isopropylalcohol.png'

Hvap = dict()
Hvap['Acetona'] = 31300
Hvap['Acetonitrilo'] = 33225
Hvap['Acido Acetico'] = 23700
Hvap['Agua'] = 40660
Hvap['Etanol'] = 38600
#Hvap['Etilenglicol'] = 53200
Hvap['Fenol'] = 58800
Hvap['isopropyl-alcohol'] = 44000

cPLiq = dict()
cPLiq['Acetona'] = 125.5
cPLiq['Acetonitrilo'] = 91.7
cPLiq['Acido Acetico'] = 123.1
cPLiq['Agua'] = 75.327
cPLiq['Etanol'] = 112.4
#cPLiq['Etilenglicol'] = 37.8
cPLiq['Fenol'] = 138.0
cPLiq['isopropyl-alcohol'] = 160.8 

cPVap = dict()
cPVap['Acetona'] = 75
cPVap['Acetonitrilo'] = 80 
cPVap['Acido Acetico'] = 63.4
cPVap['Agua'] = 37.47
cPVap['Etanol'] = 82.0
#cPVap['Etilenglicol'] = 23.0
cPVap['Fenol'] = 145.52
cPVap['isopropyl-alcohol'] = 92.4

def checkvolatility(FluidoA='Acetona', FluidoB='Agua'):
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
        
def getTeb(Fluido='Agua', P=1):
    print("Punto de ebullición del %s : %.3f °C"%(Fluido,Tsat[Fluido](P*760)))
    from IPython.display import Image
    u=Image(filename=Molecule[Fluido], width=150, height=200)
    display(u)
    
widgetTeb=interactive(getTeb, Fluido=Tsat.keys(), P=widgets.BoundedFloatText(value=1,min=0.1,max=5.0,step=0.1,description='p (atm)'))
controlsTeb = HBox(widgetTeb.children[:-1], layout = Layout(flex_flow='row wrap'))
outputTeb = widgetTeb.children[-1]

def checkTeb(FluidoA='Acetona', FluidoB='Agua'):
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

def checkRaoult(FluidoA='Acetona', FluidoB='Agua',T=100):
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

def checkTx(FluidoA='Acetona', FluidoB='Agua',Pext=1.0):
    if FluidoA==FluidoB:
        print("Seleccione dos fluidos diferentes")
        return
    P=Pext*760
    x = lambda T: (P-Psat[FluidoB](T))/(Psat[FluidoA](T)-Psat[FluidoB](T))
    y = lambda T: x(T)*Psat[FluidoA](T)/P
    plt.figure(2,figsize=(6, 4), dpi= 100)
    T = np.linspace(Tsat[FluidoA](P),Tsat[FluidoB](P))
    plt.plot([x(T) for T in T],T,color='black')
    plt.plot([y(T) for T in T],T,color='black')
    plt.xlabel("$x_{"+FluidoA+"}$")
    plt.ylabel('Temperatura $^\circ$C')
    plt.title('Diagrama Txy para {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
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

def checkFlash(FluidoA='Acetona', FluidoB='Agua',Pext=1.0,xF=0.5,TF=100.0,F=100.0):
    if FluidoA==FluidoB:
        print("Seleccione dos fluidos diferentes")
        return
    P=Pext*760
    x = lambda T: (P-Psat[FluidoB](T))/(Psat[FluidoA](T)-Psat[FluidoB](T))
    y = lambda T: x(T)*Psat[FluidoA](T)/P
    plt.figure(2,figsize=(6, 4), dpi= 100)
    T = np.linspace(Tsat[FluidoA](P),Tsat[FluidoB](P))
    plt.plot([x(T) for T in T],T,color='black')
    plt.plot([y(T) for T in T],T,color='black')
    plt.xlabel("$x_{"+FluidoA+"}$")
    plt.ylabel('Temperatura $^\circ$C')
    plt.title('Diagrama Txy para {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
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
        print("L = %.3f mol/s -- V = %.3f mol/s -- Tburbuja = %.1f °C -- Trocío = %.f °C"%(L, V, Tbub[0], Tdew[0]))
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

def checkyx(FluidoA='Acetona', FluidoB='Agua',Pext=1.0):
    if FluidoA==FluidoB:
        print("Seleccione dos fluidos diferentes")
        return
    P=Pext*760
    x = lambda T: (P-Psat[FluidoB](T))/(Psat[FluidoA](T)-Psat[FluidoB](T))
    y = lambda T: x(T)*Psat[FluidoA](T)/P
    T = np.linspace(Tsat[FluidoA](P),Tsat[FluidoB](P))
    plt.figure(figsize=(5,5),dpi=100)
    plt.plot([x(T) for T in T],[y(T) for T in T], color='black')
    plt.plot([0,1],[0,1],color='black',linestyle='--')
    plt.axis('equal')

    plt.title('Diagrama x-y para {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
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
    FluidoA='isopropyl-alcohol'
    FluidoB='Acido Acetico'
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
    plt.ylabel('Temperatura $^\circ$C')
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

    plt.suptitle('Diagramas T-xy e y-x para {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
    plt.subplots_adjust(wspace=0.2)
    plt.show()
    
widgetFeed=interactive(checkFeed,
                         TF=widgets.BoundedFloatText(value=100,min=100,max=150,step=1,description='$T_F$ (°C)'), 
                         xF=widgets.BoundedFloatText(value=0.5,min=0.1,max=0.9,step=0.01,description='$x_F$'))
controlsFeed = HBox(widgetFeed.children[:-1], layout = Layout(flex_flow='row wrap'))
outputFeed = widgetFeed.children[-1]

##### Diagrama y-x con recta de alimentación y recta de operación de la zona de rectificación

def checkFeedRectif(TF=100.0,xF=0.5,xD=0.9,Rext=3.0):
    FluidoA='isopropyl-alcohol'
    FluidoB='Acido Acetico'
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
        print('Reflujo menor que el mínimo permitido (%.4f)'%Rextmin)
    
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
    
    plt.title('Diagrama x-y para {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
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
    FluidoA='isopropyl-alcohol'
    FluidoB='Acido Acetico'
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
        print('Reflujo menor que el mínimo permitido (%.4f)'%Rextmin)
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
        
        
    plt.title('Diagrama x-y para {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
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

def checkMcCabeStep(TF=100.0,xF=0.5,xB=0.1,xD=0.9,Rext=3.0,Pasos=1):
    FluidoA='isopropyl-alcohol'
    FluidoB='Acido Acetico'
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
        print('Reflujo menor que el mínimo permitido (%.4f)'%Rextmin)
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
    
    plt.title('Diagrama x-y para {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
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

    while xP > xB and nTray < Pasos:
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
    
    if xP <= xB:
        plt.text(0.05,0.9,'{:d} etapas eq.\nRmin={:.4g}'.format(int(nTray),Rextmin),weight='bold')
    else:
        plt.text(0.05,0.9,'\nRmin={:.4g}'.format(Rextmin),weight='bold')

    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.minorticks_on()
    plt.grid(linewidth=1, which='both')
    
widgetMcCabeStep=interactive(checkMcCabeStep,
                         TF=widgets.BoundedFloatText(value=100,min=30,max=300,step=5,description='$T_F$ (°C)'), 
                         xF=widgets.BoundedFloatText(value=0.5,min=0.1,max=0.9,step=0.01,description='$x_F$'), 
                         xB=widgets.BoundedFloatText(value=0.1,min=0.01,max=0.49,step=0.01,description='$x_B$'), 
                         xD=widgets.BoundedFloatText(value=0.9,min=0.51,max=0.99,step=0.01,description='$x_D$'),                           Rext=widgets.BoundedFloatText(value=3.0,min=0.1,max=5.0,step=0.1,description='$R_{ext}$'),
                         Pasos=widgets.BoundedIntText(value=1,min=1,max=50,step=1))
controlsMcCabeStep = HBox(widgetMcCabeStep.children[:-1], layout = Layout(flex_flow='row wrap'))
outputMcCabeStep = widgetMcCabeStep.children[-1]


##### Diagrama completo McCabe-Thiele

def checkMcCabe(FluidoA='Acetona', FluidoB='Agua',Pext=1.0,TF=100.0,xF=0.5,xB=0.1,xD=0.9,Rext=3.0):
    if FluidoA==FluidoB:
        print("Seleccione dos fluidos diferentes")
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
        print("Seleccione el FluidoA como el más volátil")
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
        print('Reflujo menor que el mínimo permitido (%.4f)'%Rextmin)
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
    
    plt.title('Diagrama x-y para {:s}/{:s} a p = {:.1f} atm'.format(FluidoA,FluidoB,Pext))
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
    
    plt.text(0.05,0.9,'{:d} etapas eq.\nRmin={:.4g}'.format(int(nTray),Rextmin),weight='bold')

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
