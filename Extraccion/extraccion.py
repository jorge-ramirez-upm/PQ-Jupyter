#%matplotlib inline
import matplotlib.pyplot as plt
import ternary
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, splrep, splev
from scipy.optimize import root
from ipywidgets import interact, interactive, fixed, interact_manual, widget, widgets, Layout, HBox, VBox
import matplotlib.patches as mpl_patches

# READ THE FILE WITH THE BINODAL CURVE AS x, y, z (tab or space separated) with headers
# The Solute should be the second column
FILENAME='Water-EthylEthanoate-EthanoicAcid.dat'
data=pd.read_csv('Water-EthylEthanoate-EthanoicAcid.dat', sep='\s+|\t', engine='python')
# In case the columns need reordering, the following may be useful
columns_titles = [data.columns[1],data.columns[2],data.columns[0]]
data=data.reindex(columns=columns_titles)
# Sort data according to column 0
data.sort_values(inplace=True, by=[data.columns[0]])

# READ THE FILE WITH THE TIE LINES AS xA, yA, zA, xB, yB, zB (tab or space separated) with headers
FILENAME='Water-EthylEthanoate-EthanoicAcid.dat'
datatie=pd.read_csv('Water-EthylEthanoate-EthanoicAcid-TieLines.dat', sep='\s+|\t', engine='python')
# In case the columns need reordering, the following may be useful
columns_titles = [datatie.columns[1],datatie.columns[2],datatie.columns[0],datatie.columns[4],datatie.columns[5],datatie.columns[3]]
datatie=datatie.reindex(columns=columns_titles)

sense=["Como el reloj", "Contrario al reloj"]

"""
 _                                   
| |_ ___ _ __ _ __   __ _ _ __ _   _ 
| __/ _ \ '__| '_ \ / _` | '__| | | |
| ||  __/ |  | | | | (_| | |  | |_| |
 \__\___|_|  |_| |_|\__,_|_|   \__, |
                               |___/ 
"""

def checkternary(xA=30, xB=30, Sentido="Contrario al reloj"):

    # Two subplots: one for ternary diagram and another one for standard axes
    ax1=plt.subplot(111)

    # Set the scale of the ternary axes
    scale = 100
    figure, tax = ternary.figure(ax1, scale=scale)
    figure.set_size_inches(7, 7)

    # Draw Boundary and Gridlines
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=5)

    # Clear matplotlib axes
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    # Set the ticks in the ternary axes
    if Sentido=="Contrario al reloj":
        tax.ticks(axis='lbr', multiple=10, linewidth=1, offset=0.025)
    else:
        tax.ticks(axis='lbr', multiple=10, linewidth=1, offset=0.025, clockwise=True)

    # Set Axis labels and Title
    fontsize = 15
    offset = 0.15
    tax.right_corner_label("B", fontsize=fontsize, offset=offset, fontweight='bold')
    tax.top_corner_label("C", fontsize=fontsize, offset=offset, fontweight='bold')
    tax.left_corner_label("A", fontsize=fontsize, offset=offset, fontweight='bold')
    if Sentido=="Contrario al reloj":
        tax.left_axis_label("A", fontsize=fontsize, fontweight='bold', offset=offset)
        tax.right_axis_label("C", fontsize=fontsize, fontweight='bold',offset=offset)
        tax.bottom_axis_label("B", fontsize=fontsize, fontweight='bold',offset=offset)
    else:
        tax.left_axis_label("C", fontsize=fontsize, fontweight='bold', offset=offset)
        tax.right_axis_label("B", fontsize=fontsize, fontweight='bold',offset=offset)
        tax.bottom_axis_label("A", fontsize=fontsize, fontweight='bold',offset=offset)

    # IS THE COMPOSITION RIGHT?
    if xA + xB > 100:
        # Use legend to show some text
        handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                         lw=0, alpha=0)] * 2

        # create the corresponding number of labels (= the text you want to display)
        labels = []
        labels.append("Composición no válida")
        labels.append("Se debe cumplir: $x_A + x_B \\leq 100$")
        legend_properties = {'weight':'bold'}
        plt.legend(handles, labels, loc='best', fontsize='x-large',
                  fancybox=True, framealpha=0.5, facecolor='tab:orange',
                  handlelength=0, handletextpad=0, prop=legend_properties)

        plt.show()
        return

    # Draw the point
    xC=100-xA-xB
    P=np.array((xB, xC, xA))
    tax.scatter([P], s=25, color='m')

    if Sentido=="Como el reloj":
        P2=np.array((100-xA, 0, xA))
        tax.line(P,P2, linestyle=':', color='r')
        tax.annotate("$x_A$", (P+P2)/2, ha='left', va='bottom', color='r', fontsize=15)
        P2=np.array((xB, 100-xB, 0))
        tax.line(P,P2, linestyle=':', color='b')
        tax.annotate("$x_B$", (P+P2)/2, ha='left', va='bottom', color='b', fontsize=15)
        P2=np.array((0, xC, 100-xC))
        tax.line(P,P2, linestyle=':', color='g')
        tax.annotate("$x_C$", (P+P2)/2, ha='left', va='bottom', color='g', fontsize=15)
    else:
        P2=np.array((0, 100-xA, xA))
        tax.line(P,P2, linestyle=':', color='r')
        tax.annotate("$x_A$", (P+P2)/2, ha='left', va='bottom', color='r', fontsize=15)
        P2=np.array((xB, 0, 100-xB))
        tax.line(P,P2, linestyle=':', color='b')
        tax.annotate("$x_B$", (P+P2)/2, ha='left', va='bottom', color='b', fontsize=15)
        P2=np.array((100-xC, xC, 0))
        tax.line(P,P2, linestyle=':', color='g')
        tax.annotate("$x_C$", (P+P2)/2, ha='left', va='bottom', color='g', fontsize=15)

    # arrow_params = {'shape': 'right'}
    # plt.arrow(0, 0, 40, 40, color='r', linestyle=':', head_length=10, 
    #           **arrow_params)
    # plt.annotate(s='LOCO', xy=(20,50), xytext=(0,0), arrowprops=dict(arrowstyle='->'))

    plt.show()

widgetternary=interactive(checkternary,
                             xA=widgets.BoundedFloatText(value=30,min=0,max=100,step=1,description='$x_A$ (%)'),
                             xB=widgets.BoundedFloatText(value=30,min=0,max=100,step=1,description='$x_B$ (%)'),
                             Sentido=sense)

controlsternary = HBox(widgetternary.children[:-1]) #, layout = Layout(flex_flow='row wrap'))
outputternary = widgetternary.children[-1]

"""
 _                     
| | _____   _____ _ __ 
| |/ _ \ \ / / _ \ '__|
| |  __/\ V /  __/ |   
|_|\___| \_/ \___|_|   
"""


def checklever(F1=100, F2=100, xA1=30, xB1=30, xA2=60, xB2=60):

    # Two subplots: one for ternary diagram and another one for standard axes
    ax1=plt.subplot(111)

    # Set the scale of the ternary axes
    scale = 100
    figure, tax = ternary.figure(ax1, scale=scale)
    figure.set_size_inches(7, 7)

    # Draw Boundary and Gridlines
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=5)

    # Clear matplotlib axes
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    # Set the ticks in the ternary axes
    tax.ticks(axis='lbr', multiple=10, linewidth=1, offset=0.025)

    # Set Axis labels and Title
    fontsize = 15
    offset = 0.15
    tax.left_axis_label("A", fontsize=fontsize, fontweight='bold', offset=offset)
    tax.right_axis_label("C", fontsize=fontsize, fontweight='bold',offset=offset)
    tax.bottom_axis_label("B", fontsize=fontsize, fontweight='bold',offset=offset)

    # IS THE COMPOSITION RIGHT?
    if (xA1 + xB1 > 100) or (xA2 + xB2 > 100):
        # Use legend to show some text
        handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                         lw=0, alpha=0)] * 4

        # create the corresponding number of labels (= the text you want to display)
        labels = []
        labels.append("Composición no válida")
        labels.append("Se debe cumplir:")
        labels.append("$x_{A1} + x_{B1} \\leq 100$")
        labels.append("$x_{A2} + x_{B2} \\leq 100$")
        legend_properties = {'weight':'bold'}
        plt.legend(handles, labels, loc='best', fontsize='x-large',
                  fancybox=True, framealpha=0.5, facecolor='tab:orange',
                  handlelength=0, handletextpad=0, prop=legend_properties)

        plt.show()
        return

    # Draw the point
    xC1=100-xA1-xB1
    xC2=100-xA2-xB2
    P1=np.array((xB1, xC1, xA1))
    P2=np.array((xB2, xC2, xA2))
    tax.scatter([P1, P2], s=25, color='g')
    tax.annotate("$P_1$", P1, ha='right', color='b', fontsize=15, fontweight='bold')
    tax.annotate("$P_2$", P2, ha='right', color='b', fontsize=15, fontweight='bold')
    tax.line(P1, P2, linewidth=1.5, color='r')

    M=(F1*P1+F2*P2)/(F1+F2)
    MQ = F1+F2
    tax.scatter([M], s=25, color='m')
    tax.annotate("$M$", M, ha='right', color='b', fontsize=15, fontweight='bold')
    
    # Use legend to show some text
    handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                        lw=0, alpha=0)] * 4

    # create the corresponding number of labels (= the text you want to display)
    labels = []
    labels.append("M: %.3g mol/s"%MQ)
    labels.append("$x_{AM}$=%.3g"%M[2])
    labels.append("$x_{BM}$=%.3g"%M[0])
    labels.append("$x_{CM}$=%.3g"%M[1])
    legend_properties = {'weight':'bold'}
    plt.legend(handles, labels, loc='best', fontsize='x-large',
                fancybox=True, framealpha=0.5, facecolor='tab:green',
                handlelength=0, handletextpad=0, prop=legend_properties)
    
    
    plt.show()

widgetlever=interactive(checklever,
                        F1=widgets.BoundedFloatText(value=100,min=0,max=300,step=10,description='$F_1$ (mol/s)'),
                        F2=widgets.BoundedFloatText(value=100,min=0,max=300,step=10,description='$F_2$ (mol/s)'),
                        xA1=widgets.BoundedFloatText(value=30,min=0,max=100,step=1,description='$x_{A1}$ (%)'),
                        xB1=widgets.BoundedFloatText(value=30,min=0,max=100,step=1,description='$x_{B1}$ (%)'),
                        xA2=widgets.BoundedFloatText(value=60,min=0,max=100,step=1,description='$x_{A2}$ (%)'),
                        xB2=widgets.BoundedFloatText(value=30,min=0,max=100,step=1,description='$x_{B2}$ (%)')
)

controlslever = HBox(widgetlever.children[:-1] , layout = Layout(flex_flow='row wrap'))
outputlever = widgetlever.children[-1]

"""
 _     _                 _       _ 
| |__ (_)_ __   ___   __| | __ _| |
| '_ \| | '_ \ / _ \ / _` |/ _` | |
| |_) | | | | | (_) | (_| | (_| | |
|_.__/|_|_| |_|\___/ \__,_|\__,_|_|
"""

def checkbinodal(xA=30, xB=30):

    """ Example to test ternary diagrams and 2D graphs simultaneously
    """

    # Two subplots: one for ternary diagram and another one for standard axes
    ax1=plt.subplot(111)

    # Set the scale of the ternary axes
    scale = 100
    figure, tax = ternary.figure(ax1, scale=scale)
    figure.set_size_inches(7, 7)

    # Draw Boundary and Gridlines
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=5)

    # Clear matplotlib axes
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    # Set the ticks in the ternary axes
    tax.ticks(axis='lbr', multiple=10, linewidth=1, offset=0.025)

    # Set Axis labels and Title
    fontsize = 10
    offset = 0.15
    tax.left_axis_label(data.columns[2], fontsize=fontsize, offset=offset)
    tax.right_axis_label(data.columns[1], fontsize=fontsize, offset=offset)
    tax.bottom_axis_label(data.columns[0], fontsize=fontsize, offset=offset)

    points=data.to_numpy()*100
    # CLEAN UP THE DATA A BIT
    npts = 100
    smooth=2.0
    pointsinterp=np.zeros((npts,3))
    pointsinterp[:,0] = np.linspace(min(points[:,0]),max(points[:,0]),npts)
    sply = splrep(points[:,0],points[:,1], s=smooth)
    pointsinterp[:,1] = splev(pointsinterp[:,0], sply)
    splz = splrep(points[:,0],points[:,2], s=smooth)
    pointsinterp[:,2] = splev(pointsinterp[:,0], splz)

    # Plot the data
    #tax.scatter(points, color='m')
    tax.plot(pointsinterp, linewidth=2.0, color='k')
    xs, ys=ternary.helpers.project_sequence(pointsinterp)
    ax1.fill(xs, ys, color='tab:orange', alpha=0.5)

    # IS THE COMPOSITION RIGHT?
    if xA + xB > 100:
        # Use legend to show some text
        handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                         lw=0, alpha=0)] * 2

        # create the corresponding number of labels (= the text you want to display)
        labels = []
        labels.append("Composición no válida")
        labels.append("Se debe cumplir: $x_A + x_B \\leq 100$")
        legend_properties = {'weight':'bold'}
        plt.legend(handles, labels, loc='best', fontsize='x-large',
                  fancybox=True, framealpha=0.5, facecolor='tab:orange',
                  handlelength=0, handletextpad=0, prop=legend_properties)

        plt.show()
        return
        
    # FEED AND SOLVENT (composition and volume rate)
    F=np.array((xB, 100-xA-xB, xA))
    
    # Draw F, S, M and line
    tax.scatter([F], color='b')
    tax.annotate("M", F, ha='right', color='b', fontsize=15, fontweight='bold')

    # IS M IN THE AREA BELOW THE BINODAL?
    if F[1] > splev(F[0], sply):
        # Use legend to show some text
        handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                         lw=0, alpha=0)] * 2

        # create the corresponding number of labels (= the text you want to display)
        labels = []
        labels.append("M fuera de curva binodal")
        labels.append("Extracción imposible")
        legend_properties = {'weight':'bold'}
        plt.legend(handles, labels, loc='best', fontsize='x-large',
                  fancybox=True, framealpha=0.5, facecolor='tab:red',
                  handlelength=0, handletextpad=0, prop=legend_properties)

    else:
        # Use legend to show some text
        handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                         lw=0, alpha=0)] * 2

        # create the corresponding number of labels (= the text you want to display)
        labels = []
        labels.append("M dentro de la curva binodal")
        labels.append("Extracción posible")
        legend_properties = {'weight':'bold'}
        plt.legend(handles, labels, loc='best', fontsize='x-large',
                  fancybox=True, framealpha=0.5, facecolor='tab:green',
                  handlelength=0, handletextpad=0, prop=legend_properties)

    plt.show()

widgetbinodal=interactive(checkbinodal,
                          xA=widgets.BoundedFloatText(value=30,min=0,max=100,step=1,description='$x_{AM}$ (%)'),
                          xB=widgets.BoundedFloatText(value=30,min=0,max=100,step=1,description='$x_{BM}$ (%)'),
)
controlsbinodal = HBox(widgetbinodal.children[:-1]) #, layout = Layout(flex_flow='row wrap'))
outputbinodal = widgetbinodal.children[-1]

"""
                           _        
 _ __ ___ _ __   __ _ _ __| |_ ___  
| '__/ _ \ '_ \ / _` | '__| __/ _ \ 
| | |  __/ |_) | (_| | |  | || (_) |
|_|  \___| .__/ \__,_|_|   \__\___/ 
         |_|                        
"""

def checkreparto(xA=30, xB=30):

    """ Example to test ternary diagrams and 2D graphs simultaneously
    """

    # Two subplots: one for ternary diagram and another one for standard axes
    ax1=plt.subplot(111)

    # Set the scale of the ternary axes
    scale = 100
    figure, tax = ternary.figure(ax1, scale=scale)
    figure.set_size_inches(7, 7)

    # Draw Boundary and Gridlines
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=5)

    # Clear matplotlib axes
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    # Set the ticks in the ternary axes
    tax.ticks(axis='lbr', multiple=10, linewidth=1, offset=0.025)

    # Set Axis labels and Title
    fontsize = 10
    offset = 0.15
    tax.left_axis_label(data.columns[2], fontsize=fontsize, offset=offset)
    tax.right_axis_label(data.columns[1], fontsize=fontsize, offset=offset)
    tax.bottom_axis_label(data.columns[0], fontsize=fontsize, offset=offset)

    points=data.to_numpy()*100
    pointstie=datatie.to_numpy()*100
    pointstieA=pointstie[:,:3]
    pointstieB=pointstie[:,-3:]
    distMtie=np.zeros(len(pointstie))
    # CLEAN UP THE DATA A BIT
    npts = 100
    smooth=2.0
    pointsinterp=np.zeros((npts,3))
    pointsinterp[:,0] = np.linspace(min(points[:,0]),max(points[:,0]),npts)
    sply = splrep(points[:,0],points[:,1], s=smooth)
    pointsinterp[:,1] = splev(pointsinterp[:,0], sply)
    splz = splrep(points[:,0],points[:,2], s=smooth)
    pointsinterp[:,2] = splev(pointsinterp[:,0], splz)

    # Plot the data
    #tax.scatter(points, color='m')
    tax.plot(pointsinterp, linewidth=2.0, color='k')

    # CLEAN UP TIE LINES
    for p in pointstie:
        p[1]=splev(p[0], sply)
        p[2]=splev(p[0], splz)
        p[4]=splev(p[3], sply)
        p[5]=splev(p[3], splz)
        tax.line(p[:3], p[-3:], linewidth=1.5, color='k', linestyle="--")

    # FEED AND SOLVENT (composition and volume rate)
    M=np.array((xB, 100-xA-xB, xA))

    # Draw M
    tax.scatter([M], color='b')
    tax.annotate("M", M, ha='right', color='b', fontsize=15, fontweight='bold')

    # IS M IN THE AREA BELOW THE BINODAL?
    if M[1] > splev(M[0], sply):
        # Use legend to show some text
        handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                         lw=0, alpha=0)] * 2

        # create the corresponding number of labels (= the text you want to display)
        labels = []
        labels.append("F fuera de curva binodal")
        labels.append("Extracción imposible")
        legend_properties = {'weight':'bold'}
        plt.legend(handles, labels, loc='best', fontsize='x-large',
                  fancybox=True, framealpha=0.5, facecolor='tab:red',
                  handlelength=0, handletextpad=0, prop=legend_properties)

        plt.show()
        return

    # FIND DISTANCE FROM M TO ALL TIE LINES
    for i in range(len(pointstieA)):
        distMtie[i] = np.cross(pointstieB[i, :2]-pointstieA[i, :2],M[:2]-pointstieA[i, :2])/np.linalg.norm(pointstieB[i, :2]-pointstieA[i, :2])

    # CASES: 
    # 1) All distances are positive: Take the smallest distance and use the slope
    if np.all(distMtie>0):
        indbelow=distMtie.argmin()
        slopebelow=(pointstieB[indbelow,1]-pointstieA[indbelow,1])/(pointstieB[indbelow,0]-pointstieA[indbelow,0])
        slope=slopebelow
    # 2) All distances are negative: Take the smallest two distances (in abs) and work a slope
    elif np.all(distMtie<0):
        indabove=distMtie.argmax()
        slopeabove=(pointstieB[indabove,1]-pointstieA[indabove,1])/(pointstieB[indabove,0]-pointstieA[indabove,0])
        slope=slopeabove
    # 3) Some are positive and some negative: Take the smallest positive and the largest negative and work out a slope
    else:
        indbelow=np.where(distMtie < 0, distMtie, -np.inf).argmax()
        slopebelow=(pointstieB[indbelow,1]-pointstieA[indbelow,1])/(pointstieB[indbelow,0]-pointstieA[indbelow,0])
        indabove=np.where(distMtie > 0, distMtie, np.inf).argmin()
        slopeabove=(pointstieB[indabove,1]-pointstieA[indabove,1])/(pointstieB[indabove,0]-pointstieA[indabove,0])
        slope=(np.abs(distMtie[indbelow])*slopeabove + np.abs(distMtie[indabove])*slopebelow)/(np.abs(distMtie[indbelow])+np.abs(distMtie[indabove]))

    # Find intersection between tie line and binodal
    ysign=np.sign(M[1]+slope*(pointsinterp[:,0]-M[0])-pointsinterp[:,1])
    signchange = ((np.roll(ysign, 1) - ysign) != 0).astype(int)
    positions=np.argwhere(signchange)
    f = lambda x: M[1]+slope*(x-M[0])-splev(x, sply) 
    sol=[]
    for i in positions:
        solnow = root(f, pointsinterp[i[0],0])
        sol.append(solnow.x[0])
    A = np.array([sol[0], splev(sol[0], sply)])
    B = np.array([sol[1], splev(sol[1], sply)])

    # Draw the new tie line
    A3=[A[0], A[1], 100-A[0]-A[1]]
    B3=[B[0], B[1], 100-B[0]-B[1]]
    tax.scatter([A3, B3], color='r', marker='s', s=30)
    tax.line(A3, B3, color='r', linestyle='--')

    tax.annotate("R", A3, ha='right', va='bottom', color='r', fontsize=15, fontweight='bold')
    tax.annotate("E", B3, ha='left', va='bottom', color='r', fontsize=15, fontweight='bold')

    # Use legend to show some text
    handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                        lw=0, alpha=0)] * 2

    # create the corresponding number of labels (= the text you want to display)
    labels = []
    labels.append("R: (%5.3g, %5.3g, %5.3g)"%(A3[2], A3[0], A3[1]))
    labels.append("E: (%5.3g, %5.3g, %5.3f)"%(B3[2], B3[0], B3[1]))
    legend_properties = {'weight':'bold'}
    plt.legend(handles, labels, loc='best', fontsize='x-large',
                fancybox=True, framealpha=0.5, facecolor='tab:green',
                handlelength=0, handletextpad=0, prop=legend_properties)


    plt.show()

widgetreparto=interactive(checkreparto,
                          xA=widgets.BoundedFloatText(value=30,min=0,max=100,step=1,description='$x_A$ (%)'),
                          xB=widgets.BoundedFloatText(value=30,min=0,max=100,step=1,description='$x_B$ (%)'),
)
controlsreparto = HBox(widgetreparto.children[:-1]) #, layout = Layout(flex_flow='row wrap'))
outputreparto = widgetreparto.children[-1]


"""
           _                      _             
  _____  _| |_ _ __ __ _  ___ ___(_) ___  _ __  
 / _ \ \/ / __| '__/ _` |/ __/ __| |/ _ \| '_ \ 
|  __/>  <| |_| | | (_| | (_| (__| | (_) | | | |
 \___/_/\_\\__|_|  \__,_|\___\___|_|\___/|_| |_|
"""


def checkextraccionsimple(Fx=30, FQ=100, SQ=100):

    """ Example to test ternary diagrams and 2D graphs simultaneously
    """

    # Two subplots: one for ternary diagram and another one for standard axes
    ax1=plt.subplot(111)

    # Set the scale of the ternary axes
    scale = 100
    figure, tax = ternary.figure(ax1, scale=scale)
    figure.set_size_inches(7, 7)

    # Draw Boundary and Gridlines
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=5)

    # Clear matplotlib axes
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    # Set the ticks in the ternary axes
    tax.ticks(axis='lbr', multiple=10, linewidth=1, offset=0.025)

    # Set Axis labels and Title
    fontsize = 10
    offset = 0.15
    #tax.right_corner_label(data.columns[0], fontsize=fontsize, offset=offset)
    #tax.top_corner_label(data.columns[1], fontsize=fontsize, offset=offset)
    #tax.left_corner_label(data.columns[2], fontsize=fontsize, offset=offset)
    tax.left_axis_label(data.columns[2], fontsize=fontsize, offset=offset)
    tax.right_axis_label(data.columns[1], fontsize=fontsize, offset=offset)
    tax.bottom_axis_label(data.columns[0], fontsize=fontsize, offset=offset)

    points=data.to_numpy()*100
    pointstie=datatie.to_numpy()*100
    pointstieA=pointstie[:,:3]
    pointstieB=pointstie[:,-3:]
    distMtie=np.zeros(len(pointstie))
    # CLEAN UP THE DATA A BIT
    npts = 100
    smooth=2.0
    pointsinterp=np.zeros((npts,3))
    pointsinterp[:,0] = np.linspace(min(points[:,0]),max(points[:,0]),npts)
    sply = splrep(points[:,0],points[:,1], s=smooth)
    pointsinterp[:,1] = splev(pointsinterp[:,0], sply)
    splz = splrep(points[:,0],points[:,2], s=smooth)
    pointsinterp[:,2] = splev(pointsinterp[:,0], splz)

    # Plot the data
    #tax.scatter(points, color='m')
    tax.plot(pointsinterp, linewidth=2.0, color='k')

    # CLEAN UP TIE LINES
    for p in pointstie:
        p[1]=splev(p[0], sply)
        p[2]=splev(p[0], splz)
        p[4]=splev(p[3], sply)
        p[5]=splev(p[3], splz)
        tax.line(p[:3], p[-3:], linewidth=1.5, color='k', linestyle="--")

    # FEED AND SOLVENT (composition and volume rate)
    F=np.array((0., Fx, 100-Fx))
    S=np.array((100., 0., 0.))
    M=(FQ*F+SQ*S)/(FQ+SQ)
    MQ = FQ+SQ

    # Draw F, S, M and line
    tax.scatter([F,S,M], color='b')
    tax.annotate("F", F, ha='right', color='b', fontsize=15, fontweight='bold')
    tax.annotate("S", S, ha='left', va='top', color='b', fontsize=15, fontweight='bold')
    tax.annotate("M", M, va='bottom', color='b', fontsize=15, fontweight='bold')
    #tax.scatter([S], color='b')
    tax.line(F, S, linewidth=1.5, color='r')

    # IS M IN THE AREA BELOW THE BINODAL?
    MIsInside = True
    if M[1] > splev(M[0], sply):
        # Use legend to show some text
        handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                         lw=0, alpha=0)] * 2

        # create the corresponding number of labels (= the text you want to display)
        labels = []
        labels.append("M fuera de curva binodal")
        labels.append("Extracción imposible")
        legend_properties = {'weight':'bold'}
        plt.legend(handles, labels, loc='best', fontsize='x-large',
                  fancybox=True, framealpha=0.5, facecolor='tab:orange',
                  handlelength=0, handletextpad=0, prop=legend_properties)

        MIsInside = False
        plt.show()
        return

    # FIND DISTANCE FROM M TO ALL TIE LINES
    for i in range(len(pointstieA)):
        distMtie[i] = np.cross(pointstieB[i, :2]-pointstieA[i, :2],M[:2]-pointstieA[i, :2])/np.linalg.norm(pointstieB[i, :2]-pointstieA[i, :2])

    # CASES: 
    # 1) All distances are positive: Take the smallest distance and use the slope
    if np.all(distMtie>0):
        indbelow=distMtie.argmin()
        slopebelow=(pointstieB[indbelow,1]-pointstieA[indbelow,1])/(pointstieB[indbelow,0]-pointstieA[indbelow,0])
        slope=slopebelow
    # 2) All distances are negative: Take the smallest two distances (in abs) and work a slope
    elif np.all(distMtie<0):
        indabove=distMtie.argmax()
        slopeabove=(pointstieB[indabove,1]-pointstieA[indabove,1])/(pointstieB[indabove,0]-pointstieA[indabove,0])
        slope=slopeabove
    # 3) Some are positive and some negative: Take the smallest positive and the largest negative and work out a slope
    else:
        indbelow=np.where(distMtie < 0, distMtie, -np.inf).argmax()
        slopebelow=(pointstieB[indbelow,1]-pointstieA[indbelow,1])/(pointstieB[indbelow,0]-pointstieA[indbelow,0])
        indabove=np.where(distMtie > 0, distMtie, np.inf).argmin()
        slopeabove=(pointstieB[indabove,1]-pointstieA[indabove,1])/(pointstieB[indabove,0]-pointstieA[indabove,0])
        slope=(np.abs(distMtie[indbelow])*slopeabove + np.abs(distMtie[indabove])*slopebelow)/(np.abs(distMtie[indbelow])+np.abs(distMtie[indabove]))

    # Find intersection between tie line and binodal
    ysign=np.sign(M[1]+slope*(pointsinterp[:,0]-M[0])-pointsinterp[:,1])
    signchange = ((np.roll(ysign, 1) - ysign) != 0).astype(int)
    positions=np.argwhere(signchange)
    f = lambda x: M[1]+slope*(x-M[0])-splev(x, sply) 
    sol=[]
    for i in positions:
        solnow = root(f, pointsinterp[i[0],0])
        sol.append(solnow.x[0])
    A = np.array([sol[0], splev(sol[0], sply)])
    B = np.array([sol[1], splev(sol[1], sply)])

    AQ=MQ*np.linalg.norm(M[:2]-B)/np.linalg.norm(A-B)
    BQ=MQ-AQ

    # Draw the new tie line
    A3=[A[0], A[1], 100-A[0]-A[1]]
    B3=[B[0], B[1], 100-B[0]-B[1]]
    tax.scatter([A3, B3], color='g', marker='s', s=30)
    tax.line(A3, B3, color='g', linestyle='--')

    tax.annotate("R", A3, ha='right', va='bottom', color='r', fontsize=15, fontweight='bold')
    tax.annotate("E", B3, ha='left', va='bottom', color='r', fontsize=15, fontweight='bold')

    # Use legend to show some text
    handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                        lw=0, alpha=0)] * 3

    # create the corresponding number of labels (= the text you want to display)
    labels = []
    labels.append("M: %.3g mol/s (%5.3g, %5.3g %5.3g)"%(MQ, 100-M[0]-M[1], M[0], M[1]))
    labels.append("R: %.3g mol/s (%5.3g, %5.3g %5.3g)"%(AQ, A3[2], A3[0], A3[1]))
    labels.append("E: %.3g mol/s (%5.3g, %5.3g %5.3g)"%(BQ, B3[2], B3[0], B3[1]))
    legend_properties = {'weight':'bold'}
    plt.legend(handles, labels, loc='best', fontsize='x-large',
                fancybox=True, framealpha=0.5, facecolor='tab:green',
                handlelength=0, handletextpad=0, prop=legend_properties)


    plt.show()

widgetextraccionsimple=interactive(checkextraccionsimple,
                             Fx=widgets.BoundedFloatText(value=30,min=1,max=100,step=1,description='$x_F$ (%)'),
                             FQ=widgets.BoundedFloatText(value=150,min=10,max=300,step=10,description='F (mol/s)'),
                             SQ=widgets.BoundedFloatText(value=100,min=10,max=300,step=10,description='S (mol/s)')
)
controlsextraccionsimple = HBox(widgetextraccionsimple.children[:-1]) #, layout = Layout(flex_flow='row wrap'))
outputextraccionsimple = widgetextraccionsimple.children[-1]

"""
                 _       
 _ __ ___   __ _(_)_ __  
| '_ ` _ \ / _` | | '_ \ 
| | | | | | (_| | | | | |
|_| |_| |_|\__,_|_|_| |_|
"""

if __name__ == "__main__":
    #checkternary(xA=30, xB=20, Sentido="Como el reloj")
    #checklever(F1=100, xA1=30, xB1=30, F2=100, xA2=60, xB2=60)
    checkbinodal(xA=30, xB=30)
    #checkextraccionsimple(Fx=30, FQ=200, SQ=200)

