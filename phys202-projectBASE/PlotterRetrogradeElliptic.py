import matplotlib.pyplot as plt
from matplotlib import animation
from IPython.html.widgets import interact, fixed

import numpy as np

from moviepy.video.io.bindings import mplfig_to_npimage
import moviepy.editor as mpy

import SolverRetrogradeElliptic as Sol
from SolverRetrogradeElliptic import *


x,dx,y,dy,X,dX,Y,dY = get_solution(parabolic,Direct,Sy)

def plot():
    plt.ioff()

    fig_mpl, ax = plt.subplots(1,figsize=(8,8), facecolor='white');
    plt.xlim(-75,75);
    plt.ylim(-75,75);
    plt.sca(ax);
    ax.set_xlabel("X Position (kpc)");
    ax.set_ylabel("Y Position (kpc)");
    plt.title("Relative Motion");
    plt.legend()

    scat = ax.scatter(X[0],Y[0],color='red',label='Disrupting Galaxy, Mass S',s=10);
    scatt = ax.scatter(x,y,color='green',s=5);
    scattt = ax.scatter(0,0,color='black',label='Central Galaxy, Mass M',s=10);


    def make_gif(u):
        newX,newY = X[u*40],Y[u*40]
        scat.set_offsets(np.transpose(np.vstack([newX,newY])))
        newx=[x[k][u*40] for k in range(120)]
        newy=[y[k][u*40] for k in range(120)]
        scatt.set_offsets(np.transpose(np.vstack([newx,newy])))
        return mplfig_to_npimage(fig_mpl)
    if parabolic == 1:
        animation = mpy.VideoClip(make_gif,duration=tmax*25);
    else:
        animation = mpy.VideoClip(make_gif,duration=tmax*5);
    #animation.ipython_display(fps=50)
    animation.write_videofile("Retrograde Elliptic Passage.mp4", fps=50)



def plotcom():
    plt.ioff()

    XCOM = S*X/(S + M)
    YCOM = S*Y/(S + M)

    XS = XCOM*(M/S)
    YS = YCOM*(M/S)
    
    XM = -XCOM
    YM = -YCOM

    xnew = x - XCOM
    ynew = y - YCOM
        
    fig_mpl, ax = plt.subplots(1,figsize=(8,8), facecolor='white');
    plt.xlim(-75,75);
    plt.ylim(-75,75);
    plt.sca(ax);
    plt.title("Center of Mass Motion");
    plt.legend()
    ax.set_xlabel("X Position (kpc)");
    ax.set_ylabel("Y Position (kpc)");
    
    scat = ax.scatter(XS,YS,color='red',label = "Disrupting Galaxy, Mass S",s=10);
    scatt = ax.scatter(xnew,ynew,color='green',s=5);
    scattt = ax.scatter(XM,YM,color='black',label = "Central Galaxy, Mass M",s=10);
    
    def make_gif(u):
        newXS,newYS = XS[u*40],YS[u*40]
        newXM,newYM = XM[u*40],YM[u*40]
        scat.set_offsets(np.transpose(np.vstack([newXS,newYS])));
        newx=[xnew[k][u*40] for k in range(120)]
        newy=[ynew[k][u*40] for k in range(120)]
        scatt.set_offsets(np.transpose(np.vstack([newx,newy])));
        scattt.set_offsets(np.transpose(np.vstack([newXM,newYM])));
        return mplfig_to_npimage(fig_mpl)

    if parabolic == 1:
        animation = mpy.VideoClip(make_gif,duration=tmax*25);
    else:
        animation = mpy.VideoClip(make_gif,duration=tmax*5);
    #animation.ipython_display(fps=50)
    animation.write_videofile("Retrograde Elliptic Passage COM.mp4", fps=50)