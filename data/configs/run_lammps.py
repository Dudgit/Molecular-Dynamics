#%%
from lammps import lammps


def run(input_filename : str):
    lmp = lammps()
    lmp.file(f"{input_filename}")

if __name__ == '__main__':
    run("water520.lmp")
# %%
