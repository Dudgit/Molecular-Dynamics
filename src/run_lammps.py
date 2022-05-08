#%%
from lammps import lammps
DATAPATH = "../data/configs/"

def run(input_filename : str):
    lmp = lammps()
    lmp.file(f"{DATAPATH}{input_filename}")

if __name__ == '__main__':
    run("water520.lmp")
