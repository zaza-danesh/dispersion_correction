#!/usr/bin/env python3

#running command: python soup_analysis.py -gmx gmx_seq -make_ndx do -msd do -rotacf do -sasa do
import argparse, os, subprocess
# from pathlib import Path
import numpy as np

 

with open('../tbl_liquids_density.tex', 'w') as tbl_density:
    tbl_density.write("\\begin{tabular}{lccccccccccccc}"+"\n")
    tbl_density.write("\\hline"+"\n")
    tbl_density.write("\\multicolumn{1}{c}{FF/LJ scheme} & \\multicolumn{6}{c}{  GAFF} & \\multicolumn{6}{c}{ CGenFF} & \\multicolumn{1}{c}{ Exp.}  \\\\"+"\n")
    tbl_density.write(" &cutoff* &cutoff  & pme &pme & pme  & pme  &cutoff*  &cutoff  & pme  & pme & pme & pme & "+"\\\\"+"\n")
    tbl_density.write(" molecules & 0 & 0  & 0 &-1 & -2  & -4  & 0  &0  & 0  & -1 & -2 & -4 & "+"\\\\"+"\n")
    tbl_density.write("\\hline"+"\n")



def parseArguments():
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    rho_exp = []
    counter=-1
    with open('density_molecules.dat', 'r') as density_exp:
        for line_den in density_exp:
            rho = []
            line_den = line_den.split("|")
            compound = line_den[0].split('.sdf')[0]
            rho_exp.append(float(line_den[1]))
    with open('iupac.txt', 'r') as density:
        for line_den in density:
            counter=counter+1
            rho = []
            line_den = line_den.split("|")
            directory = line_den[0].split('.sdf')[0]
            compound=line_den[1]
            rho_calc= []
            for ff in 'GAFF', 'CGenFF':
                for disp_coeff in 6, 5, 0, 1, 2, 4:
                    if disp_coeff == 6:
                        path='../RESULTS/LiquidsBench/LJCUTOFF-DispNo/'+ff+'/'+directory+'/'
                        lj="LJ-CutOff"
                        perc=0
                    elif disp_coeff == 5:
                        path='../RESULTS/LiquidsBench/LJCUTOFF/'+ff+'/'+directory+'/'
                        lj="LJ-CutOff"
                        perc=0
                    else:
                        path='../RESULTS/LiquidsBench/LJPME/'+ff+'/'+directory+'/'
                        lj="LJ-PME"
                        perc=disp_coeff
                    with open(path+directory+'_'+str(perc)+'_density.xvg', 'r') as fn:
                        for line in fn:
                            if line.find("@") < 0 and line.find("#") < 0:
                                data = line.split()
                                rho.append(float(data[1]))
                    rho_calc.append(np.mean(rho))
            with open('../tbl_liquids_density.tex', 'a') as tbl_density:
                tbl_density.write(compound+" & " +str(int(rho_calc[0]))+ " & " +str(int(rho_calc[1]))+ " & " +str(int(rho_calc[2]))+ " & " +str(int(rho_calc[3]))+ " & " +str(int(rho_calc[4]))+ " & " +
                    str(int(rho_calc[5]))+ " & " +str(int(rho_calc[6]))+ " & " +str(int(rho_calc[7]))+ " & " +str(int(rho_calc[8]))+ " & " +
                    str(int(rho_calc[9]))+ " & " +str(int(rho_calc[10])) + " & " +str(int(rho_calc[11])) +" & " +str(int(rho_exp[counter])) +"\\\\"+"\n")
with open('../tbl_liquids_density.tex', 'a') as tbl_density:
    tbl_density.write("\\hline"+"\n")
    tbl_density.write("\\end{tabular}"+"\n")
