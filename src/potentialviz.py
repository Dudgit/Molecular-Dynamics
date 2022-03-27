import matplotlib.pyplot as plt
import numpy as np
from calculations import LJC,dFH2,LJ,LBcomb
from spce_vars import *

start = 0.4
end = 6.0
r = np.linspace(start, end, 1000)

def plotLJ(lj, ljc, r, limx, limy, name):
    plt.subplot(1,2,1)

    plt.plot(r, lj, label='LJ')
    plt.plot([start, end], [0,0], 'k:')
    plt.title('LJ potential of '+name)
    plt.xlim(limx)
    plt.ylim(limy)
    plt.legend()

    print('LJ minimumhely: ',r[np.argmin(lj)])
    print('LJ érték: ',lj[np.argmin(lj)])

    plt.subplot(1,2,2)
    plt.plot(r, ljc, label='LJC')
    plt.plot([start, end], [0,0], 'k:')
    plt.title('LJ + Coulomb potential of '+name)
    plt.legend()

    print('LJC minimumhely: ',r[np.argmin(ljc)])
    print('LJC érték: ',ljc[np.argmin(ljc)])

def plotLJCohf(LJC_fo, LJC_fh, r, name):
    #O-H~~~F
    LJC_ohf = LJC_fo[107:997] + LJC_fh[0:890]
    
    # O-H~~~F
    #H
    ff_eps_fh = 0.0
    ff_sig_fh = 0.0 #it actually does not matter so I fixed it to zero
    ff_q_f = -1.0
    LJC_ohf2 = LJC_fo[107:997] + LJC_fh[0:890] + LJC(ff_eps_fh, ff_sig_fh, ff_q_f, spce_q_h, dFH2(r[107:997]))
    
    fig = plt.figure(figsize=(12,8))
    plt.plot(r[0:890], LJC_ohf[0:890], label='neglecting the other H')
    plt.plot(r[0:890], LJC_ohf2, label='counting in the other H too')
    plt.plot([start, end], [0,0], 'k:')
    minhely1 = r[np.argmin(LJC_ohf)]
    minhely2 = r[np.argmin(LJC_ohf2)]
    minertek1 = LJC_ohf[np.argmin(LJC_ohf)]
    minertek2 = LJC_ohf2[np.argmin(LJC_ohf2)]
    plt.plot(minhely1, minertek1, 'ro', label='minimum at: '+'%1.2f' % (minhely1) +' A')
    plt.plot(minhely2, minertek2, 'mo', label='minimum at: '+'%1.2f' % (minhely2) +' A')

    plt.title('LJ+Coulomb potential felt by $F^-$ as a function of distance from, H,\n if O-H-F are collinear (ff_'+name+')', fontsize=16)
    plt.xlabel('H-F distance [A]', fontsize=14)
    plt.ylabel('potential', fontsize=14)
    plt.legend(fontsize=12)
    plt.xlim([0.7,4])
    plt.ylim([min(LJC_ohf2)-0.05,2.5])
    
    return fig

def potentials(ff_eps_na, ff_sig_na, ff_q_na, ff_eps_f, ff_sig_f, ff_q_f, r, name):
        
    LJ_ff = LJ(ff_eps_f, ff_sig_f, r)
    LJC_ff = LJC(ff_eps_f, ff_sig_f, ff_q_f, ff_q_f, r)
    
    plotLJ(LJ_ff, LJC_ff, r, [start-0.01, 2], [min(LJ_ff)-0.05, 0.005], 'F-F')

    (ff_eps_fo, ff_sig_fo) = LBcomb(ff_eps_f, spce_eps_o, ff_sig_f, spce_sig_o)
    LJ_fo = LJ(ff_eps_fo, ff_sig_fo, r)
    LJC_fo = LJC(ff_eps_fo, ff_sig_fo, ff_q_f, spce_q_o, r)
    plotLJ(LJ_fo, LJC_fo, r, [1.65,3], [min(LJ_fo)-0.05, 0.2], 'F-O')

    (ff_eps_fh, ff_sig_fh) = LBcomb(ff_eps_f, spce_eps_h, ff_sig_f, spce_sig_h)
    LJ_fh = LJ(ff_eps_fh, ff_sig_fh, r)
    LJC_fh = LJC(ff_eps_fh, ff_sig_fh, ff_q_f, spce_q_h, r)
    plotLJ(LJ_fh, LJC_fh, r, [0,6], [min(LJ_fh)-0.05, 0.2], 'F-H')
    
    fig = plotLJCohf(LJC_fo, LJC_fh, r, name)
    
    return fig