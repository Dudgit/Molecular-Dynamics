#%%
from ELTE import *

LJ_oo = LJ(spce_eps_o, spce_sig_o, r)
LJC_oo = LJC(spce_eps_o, spce_sig_o, spce_q_o, spce_q_o, r)

plotLJ(LJ_oo, LJC_oo, r, [3,6], [-0.7,0.4], 'O-O')
# %%
figJC = potentials(ff_JC_eps_na, ff_JC_sig_na, 1.0, ff_JC_eps_f, ff_JC_sig_f, -1.0, r, 'JC')
plt.savefig("../figures/LJ_test1.png")

# %%
