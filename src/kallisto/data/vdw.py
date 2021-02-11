# src/kallisto/data/vdw.py

import numpy as np

truhlar = np.zeros(shape=(86,), dtype=np.int64)

# Scaling factors to match vdw radii from Table 12
# in J. Phys. Chem. A, Vol. 113, No. 19, 2009
truhlar = [
    0.65,  # H
    1.00,  # He
    0.65,  # Li
    0.67,  # Be
    0.92,  # B
    0.89,  # C
    0.87,  # N
    0.89,  # O
    0.91,  # F
    1.00,  # Ne
    0.81,  # Na
    0.69,  # Mg
    0.76,  # Al
    0.93,  # Si
    0.84,  # P
    0.87,  # S
    0.89,  # Cl
    1.00,  # Ar
    0.88,  # K
    0.83,  # Ca
    # Sc-Zn
    # interpolated wrt Z
    # Start: Ca factor
    # End:   Ga factor
    0.82,  # Sc
    0.82,  # Ti
    0.82,  # V
    0.81,  # Cr
    0.81,  # Mn
    0.80,  # Fe
    0.80,  # Co
    0.79,  # Ni
    0.79,  # Cu
    0.79,  # Zn
    #
    0.78,  # Ga
    0.92,  # Ge
    0.84,  # As
    0.89,  # Se
    0.88,  # Br
    1.00,  # Kr
    0.96,  # Rb
    0.86,  # Sr
    # Y-Cd
    # interpolated wrt Z
    # Start: Sr factor
    # End:   In factor
    0.86,  # Y
    0.85,  # Zr
    0.84,  # Nb
    0.83,  # Mo
    0.83,  # Tc
    0.82,  # Ru
    0.81,  # Rh
    0.80,  # Pd
    0.79,  # Ag
    0.79,  # Cd
    #
    0.78,  # In
    0.90,  # Sn
    0.89,  # Sb
    0.91,  # Te
    0.90,  # I
    1.00,  # Xe
    1.04,  # Cs
    0.89,  # Ba
    0.87,  # La
    # do not scale lanthanoids
    1.00,  # Ce
    1.00,  # Pr
    1.00,  # Nd
    1.00,  # Pm
    1.00,  # Sm
    1.00,  # Eu
    1.00,  # Gd
    1.00,  # Tb
    1.00,  # Dy
    1.00,  # Ho
    1.00,  # Er
    1.00,  # Tm
    1.00,  # Yb
    1.00,  # Lu
    # Hf-Hg
    # interpolated wrt Z
    # Start: Ba factor
    # End:   Tl factor
    0.87,  # Hf
    0.86,  # Ta
    0.85,  # W
    0.84,  # Re
    0.83,  # Os
    0.82,  # Ir
    0.81,  # Pt
    0.81,  # Au
    0.80,  # Hg
    #
    0.79,  # Tl
    0.83,  # Pb
    0.88,  # Bi
    0.85,  # Po
    0.89,  # At
    1.00,  # Rn
]


rahm = np.zeros(shape=(86,), dtype=np.int64)

# Scaling factors to match vdw radii
# in Chem. Eur. J. 2017, 23, 4017
rahm = [
    0.91,  # H
    0.95,  # He
    0.80,  # Li
    0.95,  # Be
    0.98,  # B
    0.99,  # C
    1.00,  # N
    1.00,  # O
    1.00,  # F
    1.00,  # Ne
    0.79,  # Na
    0.96,  # Mg
    0.99,  # Al
    1.02,  # Si
    1.04,  # P
    1.04,  # S
    1.05,  # Cl
    1.04,  # Ar
    0.75,  # K
    0.97,  # Ca
    0.98,  # Sc
    0.97,  # Ti
    0.97,  # V
    0.92,  # Cr
    0.96,  # Mn
    0.97,  # Fe
    0.96,  # Co
    0.94,  # Ni
    0.92,  # Cu
    0.96,  # Zn
    0.97,  # Ga
    1.03,  # Ge
    1.05,  # As
    1.05,  # Se
    1.06,  # Br
    1.05,  # Kr
    0.76,  # Rb
    0.97,  # Sr
    0.99,  # Y
    1.00,  # Zr
    0.95,  # Nb
    0.97,  # Mo
    0.98,  # Tc
    0.97,  # Ru
    0.95,  # Rh
    0.88,  # Pd
    0.95,  # Ag
    1.00,  # Cd
    0.99,  # In
    1.03,  # Sn
    1.07,  # Sb
    1.07,  # Te
    1.08,  # I
    1.08,  # Xe
    0.76,  # Cs
    0.97,  # Ba
    0.98,  # La
    1.21,  # Ce
    0.97,  # Pr
    0.97,  # Nd
    0.98,  # Pm
    0.97,  # Sm
    0.98,  # Eu
    1.01,  # Gd
    0.95,  # Tb
    0.94,  # Dy
    0.97,  # Ho
    0.97,  # Er
    0.97,  # Tm
    1.00,  # Yb
    1.00,  # Lu
    1.01,  # Hf
    1.01,  # Ta
    1.02,  # W
    1.01,  # Re
    1.01,  # Os
    1.02,  # Ir
    1.01,  # Pt
    1.00,  # Au
    1.02,  # Hg
    0.97,  # Tl
    1.02,  # Pb
    1.06,  # Bi
    1.08,  # Po
    1.09,  # At
    1.09,  # Rn
]
