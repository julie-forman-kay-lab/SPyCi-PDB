#!/bin/csh

###
# steric PALES

pales -pdb 1IGD_H_s.pdb -outD ssiaA.tbl

pales -inD dc_1IGD.tab -pdb 1IGD_H_s.pdb -outD ssiaB.tbl \
      -s1 2 -a1 7 -pdbRot rot.pdb

pales -inD dc_1IGD.tab -pdb 1IGD_H_s.pdb -outD ssiaC.tbl \
      -s1 2 -a1 7 -r1 10 -rN 60 \
      -bic -wv 0.03 -dGrid 0.3 -dot 133 -digPsi 23 -rM 33.0 \
      -lcS 0.83 -rA 0.3 -H -nosurf


###
# free format PALES

pales -stPalesFree -pdbF 1IGD_H_s.pdb -outA ssiaF.tbl

pales -stPalesFree -pdbF 1IGD_H_s.pdb -outA ssiaF.tbl \
      -bic -wv 0.03 -dGrid 0.3 -dot 133 -digPsi 23 -rM 33.0 \
      -lcS 0.83 -rA 0.3 -H -nosurf


###
# best-fit

pales -bestFit -inD dc_1IGD.tab -pdb 1IGD_H_s.pdb \
      -outD svd.tbl -s1 2 -a1 7 -pdbRot rot.pdb   \
      -map 500 -outMap worldmap.txt

pales -bestFit -inD dc_1IGD.tab -pdb 1IGD_H_s.pdb \
      -outD saupePred.tbl -s1 2 -a1 7 \
      -saupe -9.2042e-05  2.3990e-04  3.8255e-04 -4.4549e-04  3.9788e-04

pales -bestFit -inD dc_1IGD.tab -pdb 1IGD_H_s.pdb \
      -outD dadrFixed.tbl -s1 2 -a1 7 \
      -fixed -da -4.135779e-04 -dr -6.901196e-05  \
      -psi -43.23 -theta 149.43 -phi 81.40

pales -bestFit -inD dc_1IGD.tab -pdb 1IGD_H_s.pdb \
      -outD dadrOnlyFixed.tbl -s1 2 -a1 7 \
      -dadr -da -4.135779e-04 -dr -6.901196e-05


###
# DC analysis

pales -anDC -outD anDC.tbl -inD1 dc_1IGD.tab -inD2 dc_1IGD.tab \
      -s1 20 -sN 40


###
# Order matrix analysis

pales -anA -outA anA.tbl \
      -inS1 -8.9631e-05  2.4300e-04  3.8479e-04 -4.4164e-04  3.9631e-04 \
      -inS2 -1.3042e-05  4.1560e-04  3.5832e-04 -4.6099e-04  4.0923e-04
