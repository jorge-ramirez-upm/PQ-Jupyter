import numpy as np
#%matplotlib inline
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
import scipy.integrate as integrate
from ipywidgets import interact, interactive, fixed, interact_manual, widget, widgets, Layout, HBox, VBox

rA = dict()
rA['Orden -1'] = lambda x, k, cA0: (k/cA0/(1-x))
rA['Orden 0'] = lambda x, k, cA0: k*np.ones_like(x)
rA['Orden 1'] = lambda x, k, cA0: (k*cA0*(1-x))
rA['Orden 2'] = lambda x, k, cA0: (k*cA0**2*(1-x)**2)
rA['Langmuir-Hinshelwood'] = lambda x, k, cA0: (k*cA0*(1-x)*(0.5+1.2*k*cA0*(1-x))**2)
rA['Cinetica X'] = lambda x, k, cA0: 1.0/(20*k*cA0*(x-0.6)**2+4*k)
rA['Cinetica Y'] = lambda x, k, cA0: 1.0/(-20*k*cA0*(x-0.4)**2+8*k)

FA0_rA = dict()
FA0_rA['Orden -1'] = lambda x, k, cA0: 1/(k/cA0/(1-x))
FA0_rA['Orden 0'] = lambda x, k, cA0: 1/k*np.ones_like(x)
FA0_rA['Orden 1'] = lambda x, k, cA0: 1/(k*cA0*(1-x))
FA0_rA['Orden 2'] = lambda x, k, cA0: 1/(k*cA0**2*(1-x)**2)
FA0_rA['Langmuir-Hinshelwood'] = lambda x, k, cA0: 1/(k*cA0*(1-x)/(0.5+1.2*k*cA0*(1-x))**2)
FA0_rA['Cinetica X'] = lambda x, k, cA0: (20*k*cA0*(x-0.6)**2+4*k)
FA0_rA['Cinetica Y'] = lambda x, k, cA0: (-20*k*cA0*(x-0.4)**2+8*k)

def rate(kinetics='Langmuir-hinshelwood', k=0.5):
    x=np.linspace(0.01,0.99,100)
    fig, ax = plt.subplots()
    r=rA[kinetics](x,k,1.0)
    plt.plot(x,r)
    plt.xlabel("$X_A$")
    plt.ylabel("$r_A$")
    plt.xlim(0,1)
    plt.show()

def FA0_rate(kinetics='Langmuir-hinshelwood', k=0.5):
    x=np.linspace(0.01,0.99,100)
    fig, ax = plt.subplots()
    r=FA0_rA[kinetics](x,k,1.0)
    plt.plot(x,r)
    plt.xlabel("$X_A$")
    plt.ylabel("$F_{A0}/r_A$")
    plt.xlim(0,1)
    plt.show()

def plot1(Reactor='CSTR', X=0.4, kinetics='Langmuir-Hinshelwood', k=0.5):
    x=np.linspace(0,0.99,100)
    fig, ax = plt.subplots()
    r=FA0_rA[kinetics](x,k,1.0)
    plt.plot(x,r)
    plt.xlabel("$X_A$")
    plt.ylabel("$F_{A0}/r_A$")
    if kinetics=='Orden -1':
        TOP=5
    else:
        TOP=20
    plt.ylim(bottom=0, top=TOP)
    plt.xlim(0,1)
    
    if (Reactor=='CSTR'):
        y=FA0_rA[kinetics](X,k,1.0)
        V=X*y
        rect = Rectangle((0, 0), X, y, linewidth=1,edgecolor='y',facecolor='y', alpha=0.5)
        ax.add_patch(rect)
        plt.annotate("$V_{CSTR}$\n%.3g L"%V,(X/2,y*1.2),horizontalalignment='center', verticalalignment='bottom')
    else:
        V=integrate.quad(lambda x: FA0_rA[kinetics](x,k,1.0), 0, X)[0]
        ix=np.linspace(0,X,100)
        iy=FA0_rA[kinetics](ix,k,1.0)
        verts = [(0, 0), *zip(ix, iy), (X, 0)]
        poly = Polygon(verts, facecolor='y', edgecolor='y', alpha=0.5)
        ax.add_patch(poly)
        plt.annotate("$V_{PFR}$\n%.3g L"%V,(X/2,iy[-1]*1.2),horizontalalignment='center', verticalalignment='bottom')

        plt.show()

def plot2(Reactor1='CSTR', Reactor2='CSTR', X1=0.4, X2=0.8, kinetics='Langmuir-Hinshelwood', k=0.5):

    x=np.linspace(0,0.99,100)
    fig, ax = plt.subplots()
    r=FA0_rA[kinetics](x,k,1.0)
    plt.plot(x,r)
    plt.xlabel("$X_A$")
    plt.ylabel("$F_{A0}/r_A$")
    if kinetics=='Orden -1':
        TOP=5
    else:
        TOP=20
    plt.ylim(bottom=0, top=TOP)
    plt.xlim(0,1)
    
    if (Reactor1=='CSTR'):
        y1=FA0_rA[kinetics](X1,k,1.0)
        V1=X1*y1
        rect = Rectangle((0, 0), X1, y1, linewidth=1,edgecolor='y',facecolor='y', alpha=0.5)
        ax.add_patch(rect)
        plt.annotate("$V_{1,CSTR}$\n%.3g L"%V1,(X1/2,y1*1.2),horizontalalignment='center', verticalalignment='bottom')
    else:
        V1=integrate.quad(lambda x: FA0_rA[kinetics](x,k,1.0), 0, X1)[0]
        ix=np.linspace(0,X1,100)
        iy=FA0_rA[kinetics](ix,k,1.0)
        verts = [(0, 0), *zip(ix, iy), (X1, 0)]
        poly = Polygon(verts, facecolor='y', edgecolor='y', alpha=0.5)
        ax.add_patch(poly)
        plt.annotate("$V_{1,PFR}$\n%.3g L"%V1,(X1/2,iy[-1]*1.2),horizontalalignment='center', verticalalignment='bottom')

    if (Reactor2=='CSTR'):
        y2=FA0_rA[kinetics](X2,k,1.0)
        V2=(X2-X1)*y2
        rect = Rectangle((X1, 0), X2-X1, y2, linewidth=1,edgecolor='c',facecolor='c', alpha=0.5)
        ax.add_patch(rect)
        plt.annotate("$V_{2,CSTR}$\n%.3g L"%V2,((X1+X2)/2,y2*1.2),horizontalalignment='center', verticalalignment='bottom')
    else:
        V2=integrate.quad(lambda x: FA0_rA[kinetics](x,k,1.0), X1, X2)[0]
        ix=np.linspace(X1,X2,100)
        iy=FA0_rA[kinetics](ix,k,1.0)
        verts = [(X1, 0), *zip(ix, iy), (X2, 0)]
        poly = Polygon(verts, facecolor='c', edgecolor='c', alpha=0.5)
        ax.add_patch(poly)
        plt.annotate("$V_{2,PFR}$\n%.3g L"%V2,((X1+X2)/2,iy[-1]*1.2),horizontalalignment='center', verticalalignment='bottom')
                
    plt.annotate("$V=V_1+V_2$=%.3g L"%(V1+V2),(0.5,TOP-0.1),horizontalalignment='center', verticalalignment='top')
        
    plt.show()

        
def plot3(Reactor1='CSTR', Reactor2='CSTR', Reactor3='CSTR', X1=0.4, X2=0.8, X3=0.9, kinetics='Langmuir-Hinshelwood', k=0.5):

    x=np.linspace(0,0.99,100)
    fig, ax = plt.subplots()
    r=FA0_rA[kinetics](x,k,1.0)
    plt.plot(x,r)
    plt.xlabel("$X_A$")
    plt.ylabel("$F_{A0}/r_A$")
    if kinetics=='Orden -1':
        TOP=5
    else:
        TOP=20
    plt.ylim(bottom=0, top=TOP)
    plt.xlim(0,1)
    
    if (Reactor1=='CSTR'):
        y1=FA0_rA[kinetics](X1,k,1.0)
        V1=X1*y1
        rect = Rectangle((0, 0), X1, y1, linewidth=1,edgecolor='y',facecolor='y', alpha=0.5)
        ax.add_patch(rect)
        plt.annotate("$V_{1,CSTR}$\n%.3g L"%V1,(X1/2,y1*1.2),horizontalalignment='center', verticalalignment='bottom')
    else:
        V1=integrate.quad(lambda x: FA0_rA[kinetics](x,k,1.0), 0, X1)[0]
        ix=np.linspace(0,X1,100)
        iy=FA0_rA[kinetics](ix,k,1.0)
        verts = [(0, 0), *zip(ix, iy), (X1, 0)]
        poly = Polygon(verts, facecolor='y', edgecolor='y', alpha=0.5)
        ax.add_patch(poly)
        plt.annotate("$V_{1,PFR}$\n%.3g L"%V1,(X1/2,iy[-1]*1.2),horizontalalignment='center', verticalalignment='bottom')

    if (Reactor2=='CSTR'):
        y2=FA0_rA[kinetics](X2,k,1.0)
        V2=(X2-X1)*y2
        rect = Rectangle((X1, 0), X2-X1, y2, linewidth=1,edgecolor='c',facecolor='c', alpha=0.5)
        ax.add_patch(rect)
        plt.annotate("$V_{2,CSTR}$\n%.3g L"%V2,((X1+X2)/2,y2*1.2),horizontalalignment='center', verticalalignment='bottom')
    else:
        V2=integrate.quad(lambda x: FA0_rA[kinetics](x,k,1.0), X1, X2)[0]
        ix=np.linspace(X1,X2,100)
        iy=FA0_rA[kinetics](ix,k,1.0)
        verts = [(X1, 0), *zip(ix, iy), (X2, 0)]
        poly = Polygon(verts, facecolor='c', edgecolor='c', alpha=0.5)
        ax.add_patch(poly)
        plt.annotate("$V_{2,PFR}$\n%.3g L"%V2,((X1+X2)/2,iy[-1]*1.2),horizontalalignment='center', verticalalignment='bottom')
      
    if (Reactor3=='CSTR'):
        y3=FA0_rA[kinetics](X3,k,1.0)
        V3=(X3-X2)*y3
        rect = Rectangle((X2, 0), X3-X2, y3, linewidth=1,edgecolor='r',facecolor='r', alpha=0.5)
        ax.add_patch(rect)
        plt.annotate("$V_{3,CSTR}$\n%.3g L"%V3,((X2+X3)/2,y3*1.2),horizontalalignment='center', verticalalignment='bottom')
    else:
        V3=integrate.quad(lambda x: FA0_rA[kinetics](x,k,1.0), X2, X3)[0]
        ix=np.linspace(X2,X3,100)
        iy=FA0_rA[kinetics](ix,k,1.0)
        verts = [(X2, 0), *zip(ix, iy), (X3, 0)]
        poly = Polygon(verts, facecolor='r', edgecolor='r', alpha=0.5)
        ax.add_patch(poly)
        plt.annotate("$V_{3,PFR}$\n%.3g L"%V3,((X2+X3)/2,iy[-1]*1.2),horizontalalignment='center', verticalalignment='bottom')
        
        
    plt.annotate("$V=V_1+V_2+V_3$=%.3g L"%(V1+V2+V3),(0.5,TOP-0.1),horizontalalignment='center', verticalalignment='top')
        
    plt.show()
    
widgetrate=interactive(rate,
          kinetics=FA0_rA.keys(),
          k=widgets.BoundedFloatText(value=0.5,min=0.1,max=5.0,step=0.1)                   
)
controlsrate = HBox(widgetrate.children[:-1], layout = Layout(flex_flow='row wrap'))
outputrate = widgetrate.children[-1]

widgetFA_rate=interactive(FA0_rate,
          kinetics=FA0_rA.keys(),
          k=widgets.BoundedFloatText(value=0.5,min=0.1,max=5.0,step=0.1)                   
)
controlsFA_rate = HBox(widgetFA_rate.children[:-1], layout = Layout(flex_flow='row wrap'))
outputFA_rate = widgetFA_rate.children[-1]

widget1=interactive(plot1,
          Reactor=['CSTR','PFR'],
          X=widgets.BoundedFloatText(value=0.2,min=0.01,max=0.99,step=0.01,description="$X_A$"),
          kinetics=FA0_rA.keys(),
          k=widgets.BoundedFloatText(value=0.5,min=0.1,max=5.0,step=0.1)                   
)
controls1 = HBox(widget1.children[:-1], layout = Layout(flex_flow='row wrap'))
output1 = widget1.children[-1]

widget2=interactive(plot2,
          Reactor1=['CSTR','PFR'],
          Reactor2=['CSTR','PFR'],
          X1=widgets.BoundedFloatText(value=0.2,min=0.01,max=0.4,step=0.01,description="$X_{A1}$"),
          X2=widgets.BoundedFloatText(value=0.6,min=0.4,max=0.85,step=0.01,description="$X_{A2}$"),
          kinetics=FA0_rA.keys(),
          k=widgets.BoundedFloatText(value=0.5,min=0.1,max=5.0,step=0.1)                   
)
controls2 = HBox(widget2.children[:-1], layout = Layout(flex_flow='row wrap'))
output2 = widget2.children[-1]

widget3=interactive(plot3,
          Reactor1=['CSTR','PFR'],
          Reactor2=['CSTR','PFR'],
          Reactor3=['CSTR','PFR'],
          X1=widgets.BoundedFloatText(value=0.2,min=0.01,max=0.4,step=0.01,description="$X_{A1}$"),
          X2=widgets.BoundedFloatText(value=0.6,min=0.4,max=0.85,step=0.01,description="$X_{A2}$"),
          X3=widgets.BoundedFloatText(value=0.9,min=0.85,max=0.99,step=0.01,description="$X_{A3}$"),
          kinetics=FA0_rA.keys(),
          k=widgets.BoundedFloatText(value=0.5,min=0.1,max=5.0,step=0.1)                   
)
controls3 = HBox(widget3.children[:-1], layout = Layout(flex_flow='row wrap'))
output3 = widget3.children[-1]