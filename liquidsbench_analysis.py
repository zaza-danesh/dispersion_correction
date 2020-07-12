#!/usr/bin/env python3

#running command: python soup_analysis.py -gmx gmx_seq -make_ndx do -msd do -rotacf do -sasa do
import argparse, os, subprocess
# from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tkinter
# from scipy import stats
import statsmodels.api as sm

# from sklearn import datasets, linear_model

# import matplotlib
# matplotlib.use('GTKAgg')  #I had to use GTKAgg for this to work, GTK threw errors
# import matplotlib.pyplot as plt  #...  and now do whatever you need...
with open('../tbl_liquids_corr.tex', 'w') as tbl:
    tbl.write("\\begin{tabular}{ccccccc}"+"\n")
    tbl.write("\\hline"+"\n")
    tbl.write("{Force Field} & {LJ} & {D.C.}& {$\Delta C_6$} & {$a$}&{$R^2$} & {MSE\\%} \\\\"+"\n")
    tbl.write("\\hline"+"\n")


def parseArguments():
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args  = parseArguments()
    COL=['b', 'g', 'r', 'c', 'm', 'k']
    MARK=['o', '^', '*', 's', 'p', '<']
    DC=['No', 'Yes', 'No', 'No', 'No', 'No', 'No', 'Yes', 'Yes','Yes','Yes','Yes']
    dc_c=-1
    cap=['a) ', 'b) ']
    for ff in 'GAFF', 'CGenFF':
        plt.figure(num=None, figsize=(6, 6), dpi=80, facecolor='w', edgecolor='k')
        counter=-1
        for lj in 'LJCUTOFF-DispNo', 'LJCUTOFF', 'LJPME': 
            dc_c=dc_c+1
            print(ff, lj)            
            if lj == 'LJPME':
                LJ='PME'
                for percentage in 0, 1, 2, 4:
                    if percentage > 0:
                        dc_c=dc_c+1
                    counter=counter+1
                    # print(percentage)
                    rho_calc, rho_exp = [], []
                    with open('density_molecules.dat', 'r') as density:
                        for line_den in density:
                            rho = []
                            line_den = line_den.split("|")
                            compound = line_den[0].split('.sdf')[0]
                            rho_exp.append(float(line_den[1]))
                            # print("exp, ", rho_exp[-1])
                            path='../RESULTS/LiquidsBench/'+lj+'/'+ff+'/'+compound+'/'
                            with open(path+compound+'_'+str(percentage)+'_density.xvg', 'r') as fn:
                                for line in fn:
                                    if line.find("@") < 0 and line.find("#") < 0:
                                        data = line.split()
                                        rho.append(float(data[1]))
                                rho_calc.append(np.mean(rho))
                                # print("calc, ", rho_calc[-1])
                    aue = np.mean(((np.array(rho_exp) - np.array(rho_calc))/np.array(rho_exp)))
                    #Linear Regression
                    model = sm.OLS(rho_calc, rho_exp)
                    results = model.fit()
                    # print(results.params)
                    with open('../tbl_liquids_corr.tex', 'a') as tbl:
                        if percentage == 0:
                            disp_coeff = str(percentage)
                            disp_coefff = disp_coeff
                        else:
                            disp_coeff = str(-percentage)+"\%"
                            disp_coefff = str(-percentage)+"%"
                        tbl.write(ff+ " & " + LJ + " & " +DC[dc_c] +" & " +disp_coeff+ " & " + str(round(results.params[0],3))+ " & " + str(round(results.rsquared,3))+ " & "+ str(round(aue*100,2))+"\\\\"+"\n")
                    plt.scatter(rho_exp, rho_calc, marker=MARK[counter], color=COL[counter])
                    x= np.arange(np.min(rho_exp), np.max(rho_exp))
                    plt.plot(x,  results.params * x, color=COL[counter], label=LJ+' '+disp_coefff)                    
            else:
                percentage=0
                counter=counter+1
                rho_calc, rho_exp = [], []
                with open('density_molecules.dat', 'r') as density:
                    for line_den in density:
                        rho = []
                        line_den = line_den.split("|")
                        compound = line_den[0].split('.sdf')[0]
                        # print(compound)
                        rho_exp.append(float(line_den[1]))
                        path='../RESULTS/LiquidsBench/'+lj+'/'+ff+'/'+compound+'/'
                        with open(path+compound+'_'+str(percentage)+'_density.xvg', 'r') as fn:
                            for line in fn:
                                if line.find("@") < 0 and line.find("#") < 0:
                                    data = line.split()
                                    rho.append(float(data[1]))
                            rho_calc.append(np.mean(rho))
                aue = np.mean(((np.array(rho_exp) - np.array(rho_calc))/np.array(rho_exp)))
                #Linear Regression
                model = sm.OLS(rho_calc, rho_exp)
                results = model.fit()
                with open('../tbl_liquids_corr.tex', 'a') as tbl:
                    if lj == 'LJCUTOFF-DispNo':
                        LJ = 'CUTOFF*'
                    else:
                        LJ = 'CUTOFF'                        
                    tbl.write(ff+ " & " +LJ+ " & " +DC[dc_c] +" & " +str(percentage)+ " & " + str(round(results.params[0],3))+ " & " + str(round(results.rsquared,3))+ " & "+ str(round(aue*100,2))+"\\\\"+"\n")
                # fig = plt.figure()
                plt.scatter(rho_exp, rho_calc, marker=MARK[counter], color=COL[counter])
                x= np.arange(np.min(rho_exp), np.max(rho_exp))
                # print(LJ, results.params * x)
                # print(counter)
                plt.plot(x, results.params * x, color=COL[counter], label=LJ+' '+str(percentage))
                # size = fig.get_size_inches()*fig.dpi # size in pixels
                # print(size)
        if ff == 'GAFF':
            plt.text(1700, 700, '(a)', fontsize=20)
        else:
            plt.text(1700, 700, '(b)', fontsize=20)
        plt.legend(loc=2, fontsize=15)
        plt.xlim([550, 1850])
        plt.xticks([700, 900, 1100, 1300, 1500, 1700],fontsize=15)
        plt.xlabel('xlabel', fontsize=20)
        plt.ylim([550, 1850])
        # plt.ylim([500, 1900])
        plt.yticks([700, 900, 1100, 1300, 1500, 1700], fontsize=15)                
        plt.xlabel('density (exp.)', fontsize=20)
        if ff == 'GAFF':
            plt.ylabel('density (calc.)', fontsize=20)
        elif ff == 'CGenFF':
            plt.ylabel(' ', fontsize=20)
        

        # plt.axis('scaled')
        # plt.axes().set_aspect('equal', adjustable='datalim')
        plt.savefig('../fig_fit_'+ff+'.pdf', bbox_inches='tight')
        plt.close()

with open('../tbl_liquids_corr.tex', 'a') as tbl:
    tbl.write("\\hline"+"\n")
    tbl.write("\\end{tabular}"+"\n")


#     for ff in 'GAFF', 'CGenFF':
#         for counter in 6, 5, 0, 1, 2, 4:
#             rho_calc, rho_exp = [], []
#             with open('density_molecules.dat', 'r') as density:
#                 for line_den in density:
#                     rho = []
#                     line_den = line_den.split("|")
#                     compound = line_den[0].split('.sdf')[0]
#                     rho_exp.append(float(line_den[1]))
#                     if counter == 5:
#                         path='../RESULTS/LiquidsBench/LJCUTOFF/'+ff+'/'+compound+'/'
#                         lj="LJ-CutOff"
#                         perc=0
#                     elif counter == 6:
#                         path='../RESULTS/LiquidsBench/LJCUTOFF-DispNo/'+ff+'/'+compound+'/'
#                         lj="LJ-CutOff*"
#                         perc=0
#                     else:
#                         path='../RESULTS/LiquidsBench/LJPME/'+ff+'/'+compound+'/'
#                         lj="LJ-PME"
#                         perc=counter
#                     with open(path+compound+'_'+str(perc)+'_density.xvg', 'r') as fn:
#                         for line in fn:
#                             if line.find("@") < 0 and line.find("#") < 0:
#                                 data = line.split()
#                                 rho.append(float(data[1]))
#                     rho_calc.append(np.mean(rho))
#                 aue = np.mean(((np.array(rho_exp) - np.array(rho_calc))/np.array(rho_exp)))
#                 #Linear Regression
#                 model = sm.OLS(rho_calc, rho_exp)
#                 results = model.fit()
#                 print(results.params[0])
#                 with open('../tbl_liquids_corr.tex', 'a') as tbl:
#                     tbl.write(lj+' '+ff+'-'+str(perc)+ " & " + str(round(results.params[0],3))+ " & " + str(round(results.rsquared,3))+ " & "+ str(round(aue*100,2))+"\\\\"+"\n")

#                 plt.scatter(rho_exp, rho_calc, marker=MARK[counter], color=COL[counter])
#                 x= np.arange(np.min(rho_exp), np.max(rho_exp))
#                 plt.plot(x,  results.params * x, color=COL[counter], label=lj+' '+str(perc)+'%')
#                 plt.legend()
#                 plt.xlim( [np.min(rho_exp), np.max(rho_exp)])
#                 plt.xticks(fontsize=15)
#                 plt.xlabel('xlabel', fontsize=20)
#                 plt.ylim( [ min(min(rho_exp), min(rho_calc)) , max(max(rho_exp), max(rho_calc))] )                
#                 plt.yticks(fontsize=15)                
#                 plt.ylabel('ylabel',fontsize=20)
#                 plt.xlabel('density (exp.)', fontsize=20)
#                 plt.ylabel('density (calc.)', fontsize=20)
#                 plt.axis('scaled')
#         plt.savefig('../fig_fit_'+ff+'.pdf', bbox_inches='tight')
#         plt.close()



# with open('../tbl_liquids_corr.tex', 'a') as tbl:
#     tbl.write("\\hline"+"\n")
#     tbl.write("\\end{tabular}"+"\n")
