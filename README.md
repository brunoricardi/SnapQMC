# SnapQMC
1. snapshot.f90 creates snapshots of a triangular lattice forming a rhombus (U=0 and T=0).
2. snapshot_hex.f90 creates snapshots of a triangular lattice forming an hexagon (U=0 and T=0).
The outputs of each code are snapshot(_hex).out and snapcg(_hex).in. The former contains all the snapshots separated by it's probability and the final results that contains the average density and average total spin of each site. The latter contains only the snapshots and is used as an input to (h/d)cluster(_hex).f90.
3. H0.m is written on MatLab/Octave language and creates the input file psi.in to be used in snapshot.f90.
4. H_hex.m is written on MatLab/Octave language and creates the input file psi_hex.in to be used in snapshot_hex.f90.
5. (h/d)cluster(_hex).f90 reads snapcfg(_hex).in and computes prints in chist.out the frequency in which each cluster sizer appears among the snapshots; d = clusters made of doublons; h = clusters made of holes; empty = clusters made of both doublons and holes.

The code compile well with ifort. There are some warnings if compiled with gfortran.


