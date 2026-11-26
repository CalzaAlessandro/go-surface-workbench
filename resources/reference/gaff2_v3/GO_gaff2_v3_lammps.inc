# Standalone LAMMPS reference settings & 18-type LJ table for GO_gaff2_v3
# Intended for use alongside GO Workbench exported .data files (atom_style full).
units           real
atom_style      full
boundary        p p p
bond_style      harmonic
angle_style     harmonic
dihedral_style  fourier
improper_style  cvff
special_bonds   amber
pair_modify     mix arithmetic
# Pair Coeffs table (units real: epsilon [kcal/mol], sigma [Angstrom]):
#   type=ga   elem=C mass= 12.010 eps=  0.0988 sigma=  3.315212 (Rmin/2= 1.8606)
#   type=gb   elem=C mass= 12.010 eps=  0.0988 sigma=  3.315212 (Rmin/2= 1.8606)
#   type=gh   elem=H mass=  1.008 eps=  0.0161 sigma=  2.625479 (Rmin/2= 1.4735)
#   type=g6   elem=C mass= 12.010 eps=  0.1078 sigma=  3.397710 (Rmin/2= 1.9069)
#   type=gm   elem=O mass= 16.000 eps=  0.0930 sigma=  3.242871 (Rmin/2= 1.8200)
#   type=gw   elem=H mass=  1.008 eps=  0.0047 sigma=  0.537925 (Rmin/2= 0.3019)
#   type=gx   elem=C mass= 12.010 eps=  0.1078 sigma=  3.397710 (Rmin/2= 1.9069)
#   type=gy   elem=C mass= 12.010 eps=  0.1078 sigma=  3.397710 (Rmin/2= 1.9069)
#   type=gz   elem=C mass= 12.010 eps=  0.1078 sigma=  3.397710 (Rmin/2= 1.9069)
#   type=ge   elem=C mass= 12.010 eps=  0.1078 sigma=  3.397710 (Rmin/2= 1.9069)
#   type=gp   elem=O mass= 16.000 eps=  0.0726 sigma=  3.156098 (Rmin/2= 1.7713)
#   type=g1   elem=H mass=  1.008 eps=  0.0208 sigma=  2.421997 (Rmin/2= 1.3593)
#   type=gk   elem=C mass= 12.010 eps=  0.0988 sigma=  3.315212 (Rmin/2= 1.8606)
#   type=gq   elem=O mass= 16.000 eps=  0.1463 sigma=  3.048121 (Rmin/2= 1.7107)
#   type=gc   elem=C mass= 12.010 eps=  0.0988 sigma=  3.315212 (Rmin/2= 1.8606)
#   type=go   elem=O mass= 16.000 eps=  0.1463 sigma=  3.048121 (Rmin/2= 1.7107)
#   type=gd   elem=C mass= 12.010 eps=  0.0988 sigma=  3.315212 (Rmin/2= 1.8606)
#   type=gu   elem=O mass= 16.000 eps=  0.1463 sigma=  3.048121 (Rmin/2= 1.7107)
