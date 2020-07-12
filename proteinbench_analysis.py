#!/usr/bin/env python3
import argparse, os, subprocess
import numpy as np
import matplotlib.pyplot as plt
import tkinter
from matplotlib.collections import LineCollection
import matplotlib.image as image

def split_list(a_list):
    half = int(len(a_list)/2)    
    return a_list[:half], a_list[half:]

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-traj_list", "--traj_list", help="",   type=str,  default=None)
    parser.add_argument("-energetics_plot", "--energetics_plot", help="",   type=str,  default=None)
    parser.add_argument("-mindist_plot", "--mindist_plot", help="",   type=str,  default=None)
    parser.add_argument("-mindist_time", "--mindist_time", help="",   type=str,  default=None)
    parser.add_argument("-sasa_plot", "--sasa_plot", help="",   type=str,  default=None)
    parser.add_argument("-gyrate_plot", "--gyrate_plot", help="",   type=str,  default=None)
    parser.add_argument("-rmsdist_plot", "--rmsdist_plot", help="",   type=str,  default=None)
    parser.add_argument("-rmsdist_time", "--rmsdist_time", help="",   type=str,  default=None)
    parser.add_argument("-rmsf_plot", "--rmsf_plot", help="",   type=str,  default=None)
    parser.add_argument("-rdf_plot", "--rdf_plot", help="",   type=str,  default=None)
    parser.add_argument("-hbond_plot", "--hbond_plot", help="",   type=str,  default=None)
    parser.add_argument("-toc", "--toc", help="",   type=str,  default=None)

    args = parser.parse_args()
    return args

LBL_FONTSIZE=22
TCK_FONTSIZE=18
MK_SIZE=12
LGN_FONTSIZE=13
if __name__ == '__main__':
    args  = parseArguments()
    traj_directory, md_name = [], []
    MARK=['o', '^',  's', 'p', '<', '*']
    COL=['b', 'c', 'r', 'k', 'g',  'c', 'r', 'k', 'g', 'c', 'r', 'k', 'g']
    linestyle = ['-','-','-','-','-', '--','--','--','--', ':', ':', ':', ':', ':']
    LAB=['1', '2', '4', '8']



    ff_labels = []
    if args.traj_list == None:
        with open('traj_list.txt') as fn:
            for line in fn:
                line_split = line.split("/")
                traj_directory.append(line.strip())
        with open('traj_list_label.txt') as fn:
            for line in fn:
                line_split = line.split("/")
                ff_labels.append(line_split[0])
    else:
        with open(args.traj_list+'.txt') as fn:
            for line in fn:
                line_split = line.split("/")
                traj_directory.append(line.strip())
        with open(args.traj_list+'_label.txt') as fn:
            for line in fn:
                line_split = line.split("/")
                ff_labels.append(line_split[0])

    # for i in range(0, len(traj_directory)):
    #     ff_labels.append(str(traj_directory[i]))


    # path=''
    path='../RESULTS/ProteinBench/'
    if args.energetics_plot:
        if args.energetics_plot == "density":
            unit = "Density ($kg/m^3$)"
            # print(unit)
        elif args.energetics_plot == "potential":
            unit = "Potential($kJ/mol$)"        
            # print(unit)
        elif args.energetics_plot == "lj":
            unit = "Lennard Jones E.($kJ/mol$)"
            # print(unit)
        else:
            print('error')
        x = list(range(1, len(traj_directory)+1))
        c=-1
        for j in 1, 2, 4, 8:
            c=c+1
            # print('j', j)
            md_name ="ub"+str(j)
            ff_y = []
            ff_std = []
            for i in range(0, len(traj_directory)):
                # print(traj_directory[i])
                y_rep = []
                for rep in 1, 2, 3:
                    # print('rep', rep)
                    y=[]
                    with open(path+traj_directory[i]+md_name+'_'+args.energetics_plot+'.'+str(rep)+'.xvg') as fn:
                        for line in fn:
                            value = 0
                            if line.find("@") < 0 and line.find("#") < 0:
                                data = line.split()
                                y.append(float(data[1]))
                    y_rep.append(np.mean(y))
                ff_y.append(np.mean(y_rep))    
                ff_std.append(np.std(y_rep))
            plt.errorbar(x, ff_y, ff_std, fmt=MARK[c], color=COL[c], markersize=MK_SIZE, label=LAB[c])
            plt.xlim([0, len(traj_directory)+1])
            plt.xticks(x, ff_labels, rotation=80,fontsize=TCK_FONTSIZE)
            plt.yticks(fontsize=TCK_FONTSIZE)
            plt.ylabel(unit, fontsize=LBL_FONTSIZE)
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0., fontsize=TCK_FONTSIZE)

        plt.savefig('../fig_protein_'+args.energetics_plot+'.pdf', bbox_inches='tight')
        plt.close()

    if args.mindist_time:
        if args.mindist_time == "mindist":
            unit = "Min distance ($nm$)"
        elif args.mindist_plot == "numcont":
            unit = "Number of contacts (dist < 0.6 nm )"
        x = list(range(1, len(traj_directory)+1))
        j=2
        md_name ="ub"+str(j)                
        # fig = plt.figure(figsize=(8.27, 11.69))
        T = np.arange(0,310, 50)
        # print(T)
        for i in range(0, len(traj_directory)):
            fig = plt.figure(figsize=(12,0.8))
            Y = []
            for rep in 1, 2, 3:                
                y=[]
                t=[]
                with open(path+traj_directory[i]+md_name+'_'+args.mindist_time+'.'+str(rep)+'.xvg') as fn:
                    for line in fn:
                        if line.find("@") < 0 and line.find("#") < 0:
                            data = line.split()
                            if float(data[0]) > 00000:
                                y.append(float(data[1]))
                                t.append(float(data[0]))
                    Y.append(y)                               
                ax = fig.add_subplot(111)
                ax.plot(t, y)
                plt.xlim([0, 300000])
                # plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1)
            if i != 12:
                plt.xticks([],fontsize=TCK_FONTSIZE)
            else:
                plt.xticks( fontsize=TCK_FONTSIZE)
                ax.set_xticklabels([0, 50, 100, 150, 200, 250, 300])
                plt.xlabel('Time (ns)', fontsize=LBL_FONTSIZE)
            plt.ylim([0, 3])
            plt.yticks([0, 1, 2, 3], fontsize=12)
            # if i == 6:
            #     plt.ylabel(unit, fontsize=LBL_FONTSIZE)
            plt.savefig('../fig_mindist_time_'+str(i)+'.pdf', bbox_inches='tight')
            plt.close()

    if args.mindist_plot:
        if args.mindist_plot == "mindist":
            unit = "Min distance ($nm$)"
        elif args.mindist_plot == "numcont":
            unit = "Number of contacts (dist < 0.6 nm )"
        x = list(range(1, len(traj_directory)+1))
        j=2
        md_name ="ub"+str(j)                
        ff_y = []
        ff_std = []
        for i in range(0, len(traj_directory)):
            y_rep = []
            for rep in 1, 2, 3:                
                y=[]
                with open(path+traj_directory[i]+md_name+'_'+args.mindist_plot+'.'+str(rep)+'.xvg') as fn:
                    for line in fn:
                        if line.find("@") < 0 and line.find("#") < 0:
                            data = line.split()
                            if float(data[0]) > 100000:
                                y.append(float(data[1]))
                    y_rep.append(np.mean(y))
            
            ff_y.append(np.mean(y_rep))
            ff_std.append(np.std(y_rep))
        # plt.plot(x, all_g, label=j)
        # plt.errorbar(x, ff_y, ff_std, fmt='o', label=j)
        # plt.legend()
        plt.errorbar(x, ff_y, ff_std, fmt='o', markersize=MK_SIZE, color='b')
        plt.xlim([0, len(traj_directory)+1])
        plt.ylim([0, 2])
        plt.xticks(x, ff_labels, rotation=80, fontsize=TCK_FONTSIZE)
        plt.ylabel(unit, fontsize=LBL_FONTSIZE)
        plt.yticks([0, 1, 2], fontsize=TCK_FONTSIZE)
        plt.savefig('../fig_protein_'+args.mindist_plot+'.pdf', bbox_inches='tight')
        plt.close()

    if args.gyrate_plot:
        if args.gyrate_plot == "gyrate":
            unit = "Gyration radius ($nm$)"
        elif args.gyrate_plot == "":
            unit = ""
        x = list(range(1, len(traj_directory)+1))
        j=1
        md_name ="ub"+str(j)                
        ff_y = []
        ff_std = []
        for i in range(0, len(traj_directory)):
            y_rep = []
            for rep in 1, 2, 3:                
                y=[]
                with open(path+traj_directory[i]+md_name+'_'+args.gyrate_plot+'.'+str(rep)+'.xvg') as fn:
                    for line in fn:
                        if line.find("@") < 0 and line.find("#") < 0:
                            data = line.split()
                            if float(data[0]) > 100000:
                                y.append(float(data[1]))
                    y_rep.append(np.mean(y))
            
            ff_y.append(np.mean(y_rep))
            ff_std.append(np.std(y_rep))
        # plt.errorbar(x, ff_y, ff_std, fmt='o', label=j)
        # plt.legend()
        plt.errorbar(x, ff_y, ff_std, fmt='o', markersize=MK_SIZE, color='b')
        plt.xlim([0, len(traj_directory)+1])
        plt.ylim([1.165, 1.2])
        plt.xticks(x, ff_labels, rotation=80,fontsize=TCK_FONTSIZE)
        plt.yticks([1.17, 1.18, 1.19, 1.2], fontsize=TCK_FONTSIZE)
        plt.ylabel(unit, fontsize=LBL_FONTSIZE)
        plt.savefig('../fig_protein_'+args.gyrate_plot+'.pdf', bbox_inches='tight')
        plt.close()

    if args.rmsdist_plot:
        if args.rmsdist_plot == "rmsdist":
            unit = "RMSD ($nm$)"
        elif args.rmsdist_plot == "rmsdist_fromgro":
            unit = "RMSD ($nm$)"
        x = list(range(1, len(traj_directory)+1))
        j=1
        md_name ="ub"+str(j)                
        ff_y = []
        ff_std = []
        for i in range(0, len(traj_directory)):
            y_rep = []
            for rep in 1, 2, 3:                
                y=[]
                with open(path+traj_directory[i]+md_name+'_'+args.rmsdist_plot+'.'+str(rep)+'.xvg') as fn:
                    for line in fn:
                        if line.find("@") < 0 and line.find("#") < 0:
                            data = line.split()
                            if float(data[0]) > 00000:
                                y.append(float(data[1]))
                    y_rep.append(np.mean(y))
            
            ff_y.append(np.mean(y_rep))
            ff_std.append(np.std(y_rep))
        # plt.errorbar(x, ff_y, ff_std, fmt='o', label=j)
        # plt.legend()
        plt.errorbar(x, ff_y, ff_std, fmt='o', markersize=MK_SIZE, color='b')
        plt.xlim([0, len(traj_directory)+1])
        plt.ylim([0.06, 0.18])
        plt.xticks(x, ff_labels, rotation=80,fontsize=TCK_FONTSIZE)
        plt.yticks(fontsize=TCK_FONTSIZE)
        plt.ylabel(unit, fontsize=LBL_FONTSIZE)
        plt.savefig('../fig_protein_'+args.rmsdist_plot+'.pdf', bbox_inches='tight')
        plt.close()

    if args.rmsdist_time:
        if args.rmsdist_time == "rmsdist":
            unit = "RMSD ($nm$)"
        x = list(range(1, len(traj_directory)+1))
        j=1
        md_name ="ub"+str(j)                
        fig = plt.figure(figsize=(8.27, 11.69))
        T = np.arange(0,310, 50)
        # print(T)
        for i in range(0, len(traj_directory)):
            Y = []
            for rep in 1, 2, 3:                
                y=[]
                t=[]
                with open(path+traj_directory[i]+md_name+'_'+args.rmsdist_time+'.'+str(rep)+'.xvg') as fn:
                    for line in fn:
                        if line.find("@") < 0 and line.find("#") < 0:
                            data = line.split()
                            if float(data[0]) > 00000:
                                y.append(float(data[1]))
                                t.append(float(data[0]))
                    Y.append(y)               
                ax = fig.add_subplot(13, 1, i+1)
                # print(np.std(Y,axis=0))
                # plt.errorbar(t, np.mean(Y,axis=0), np.std(Y,axis=0))
                plt.plot(t, y)
                plt.title(ff_labels[i], fontsize=15)
                plt.ylim([0.05, 0.35])
                plt.xlim([0, 300000])
                plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=None, hspace=1)
            if i != 12:
                plt.xticks([],fontsize=TCK_FONTSIZE)
            else:
                plt.xticks(fontsize=TCK_FONTSIZE)
                ax.set_xticklabels([0, 50, 100, 150, 200, 250, 300])
                plt.xlabel('Time (ns)', fontsize=LBL_FONTSIZE)
            plt.yticks([0.05, 0.35], fontsize=12)
            if i == 6:
                plt.ylabel(unit, fontsize=LBL_FONTSIZE)
        plt.savefig('../fig_rmsdist_time.pdf')
        plt.close()

    if args.hbond_plot:
        if args.hbond_plot == "hb_num_pp":
            unit = "RMSD ($nm$)"
        elif args.hbond_plot == "":
            unit = ""
        x = list(range(1, len(traj_directory)+1))
        j=2
        md_name ="ub"+str(j)                
        ff_y = []
        ff_std = []
        for i in range(0, len(traj_directory)):
            y_rep = []
            for rep in 1, 2, 3:                
                y=[]
                with open(path+traj_directory[i]+md_name+'_'+args.hbond_plot+'.'+str(rep)+'.xvg') as fn:
                    for line in fn:
                        if line.find("@") < 0 and line.find("#") < 0:
                            data = line.split()
                            if float(data[0]) > 100000:
                                y.append(float(data[1]))
                    y_rep.append(np.mean(y))
            ff_y.append(np.mean(y_rep))
            ff_std.append(np.std(y_rep))
        # plt.errorbar(x, ff_y, ff_std, fmt='o', label=j)
        # plt.legend()
        plt.errorbar(x, ff_y, ff_std, fmt='o')
        plt.xlim([0, len(traj_directory)+1])
        plt.xticks(x, ff_labels, rotation=80,fontsize=TCK_FONTSIZE)
        plt.yticks(fontsize=12)
        plt.ylabel(unit, fontsize=LBL_FONTSIZE)
        plt.savefig('../fig_protein_'+args.hbond_plot+'.pdf', bbox_inches='tight')
        plt.close()

    if args.rmsf_plot:
        MK_SIZE=3
        if args.rmsf_plot == "rmsf":
            unit = "RMSF ($nm$)"
        x = list(range(1, len(traj_directory)+1))
        for i in range(0, len(traj_directory)):
            ff_y = []
            ff_std = []
            j = 1
            md_name ="ub"+str(j)                
            y_rep = []
            for rep in 1, 2, 3:                
                y=[]
                with open(path+traj_directory[i]+md_name+'_'+args.rmsf_plot+'.'+str(rep)+'.xvg') as fn:
                    for line in fn:
                        if line.find("@") < 0 and line.find("#") < 0:
                            data = line.split()
                            y.append(float(data[1]))            
                y_rep.append(y)
            plt.errorbar(list(range(1, 77)), np.mean(y_rep, axis=0), np.std(y_rep, axis=0), fmt=MARK[i], markersize=MK_SIZE, color=COL[i], label=ff_labels[i])
            # markerfacecolor='white', 
            plt.legend(loc=2, fontsize=LBL_FONTSIZE)
            plt.legend(markerscale=2)
            # lgnd = plt.legend(loc=2, fontsize=LBL_FONTSIZE)
            # for handle in lgnd.legendHandles:
            #     handle._sizes([10.0])
            # lgnd.legendHandles[0]._sizes = [30] 
            # lgnd.legendHandles[1]._sizes = [30] 
            plt.xlabel('Residue', fontsize=LBL_FONTSIZE)
            plt.ylabel(unit, fontsize=LBL_FONTSIZE)
            # plt.ylim([0.0, 0.5])
            plt.yticks(fontsize=TCK_FONTSIZE)
            plt.xticks(fontsize=TCK_FONTSIZE)
        plt.savefig('../fig_protein_'+args.rmsf_plot+'.pdf', bbox_inches='tight')
        plt.close()


    if args.rdf_plot:
        #caption=['a', 'd', 'g', 'b', 'e', 'h', 'c', 'f', 'i']
        caption=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']
        trans = [1]
        unit = "g(r)"
        if args.rdf_plot == "rdf0":
            fig_name= "rdf"
            file_ext="rdf_b100_e278ns"
        elif args.rdf_plot == "rdf1":
            fig_name= "rdf_first_half"
            file_ext="rdf_100ns"
        elif args.rdf_plot == "rdf2":
            fig_name= "rdf_second_half"
            file_ext="rdf_b200ns"
        else:
            print('Wrong argument for rdf_plot.')
        # elif args.rdf_plot == "rdf-cn":
        #     unit = "Number"
        list_ff=['amws','amcf', 'ampme', 'chpme']
        list_label=['0%', '4%']
        # list_label=['0%', '1%', '2%', '4%']
        c=0
        plt.figure(figsize=(12, 12))
        for j in  2, 4, 8:
            for I in range(4):
                traj_directory=[]
                with open(list_ff[I]+'.txt') as fn:
                # print(c)
                    for line in fn:
                        line_split = line.split("/")
                        traj_directory.append(line.strip())
                    c=c+1
                for i in range(0, len(traj_directory)):
                    # print(traj_directory[i])
                    md_name ="ub"+str(j)                
                    y_rep = []
                    all_rep = []
                    for rep in 1, 2, 3:                
                        y=[]
                        r=[]
                        with open(path+traj_directory[i]+md_name+'_'+file_ext+'.'+str(rep)+'.xvg') as fn:
                            for line in fn:
                                if line.find("@") < 0 and line.find("#") < 0:
                                    data = line.split()
                                    y.append(float(data[1]))
                                    r.append(float(data[0]))
                        y_rep.append(y)
                    m = min(len(y_rep[0]), len(y_rep[1]), len(y_rep[2]))
                    plt.subplot(3, 4, c)
                    # plt.text()
                    if I == 0:
                        ii=i
                    else:
                        ii=i+1
                    plt.plot(r[0:m], np.mean((y_rep[0][0:m],y_rep[1][0:m], y_rep[2][0:m]), axis=0),color=COL[ii], label=list_label[i])
                    if (j == 2 and I == 1):
                        plt.legend(loc=2)
                    plt.xlim([2, 4])
                    plt.xticks([2, 2.5, 3, 3.5, 4],fontsize=12)
                    plt.yticks(fontsize=12)
                if j==2:
                    plt.ylim([0, 12.5])
                    if args.rdf_plot == "rdf-cn":
                        plt.ylim([0, 2.5])

                    # plt.text(3.75, 11, str(j)+list_ff[I], fontsize=18)            
                    plt.text(3.75, 11, caption[c-1], fontsize=18)            
                elif j==4:
                    plt.ylim([0, 10.5])
                    if args.rdf_plot == "rdf-cn":
                        plt.ylim([0, 4.5])
                    # plt.text(3.75, 9, str(j)+list_ff[I], fontsize=18)            
                    plt.text(3.75, 9.1, caption[c-1], fontsize=18)            
                elif j==8:
                    plt.ylim([0, 4.5])
                    plt.yticks([0, 1, 2, 3, 4])
                    # plt.text(3.75, 4, str(j)+list_ff[I], fontsize=18)            
                    plt.text(3.75, 4, caption[c-1], fontsize=18)            

        plt.text(4.3, 14, '2 proteins', fontsize=20, rotation=90)            
        plt.text(4.3, 8.5, '4 proteins', fontsize=20, rotation=90)            
        plt.text(4.3, 3, '8 proteins', fontsize=20, rotation=90)            

        plt.text(-4.8, 15.5, 'am99-ws', fontsize=20)            
        plt.text(-2.4, 15.5, 'am99-cf', fontsize=20)            
        plt.text(-0.1, 15.5, 'am99-pme', fontsize=20)            
        plt.text(2.3, 15.5, 'ch27-pme', fontsize=20)            
        
        plt.text(-1, -1.2, 'r ($nm$)', fontsize=20)            
        plt.text(-6, 8, 'g(r)', fontsize=20, rotation=90)            
        plt.savefig('../fig_protein_'+fig_name+'.pdf',  bbox_inches='tight')
        # plt.savefig('../fig_protein_'+args.rdf_plot+'.pdf',  bbox_inches='tight')
    

    if args.toc:
        y_rep = []
        for rep in 1, 2, 3:
            y, r =[], []
            with open(path+'ffam-ws/ub8_rdf.'+str(rep)+'.xvg') as fn:
                for line in fn:
                    if line.find("@") < 0 and line.find("#") < 0:
                        data = line.split()
                        y.append(float(data[1]))
                        r.append(float(data[0]))
                y_rep.append(y)
        m = min(len(y_rep[0]), len(y_rep[1]), len(y_rep[2]))
        x = np.array(r[1100:m])
        Y = np.mean((y_rep[0][1100:m],y_rep[1][1100:m], y_rep[2][1100:m]), axis=0)        
        points = np.array([x, Y]).T.reshape(-1, 1, 2)
        segments = np.hstack([points[:-1], points[1:]])
        im = image.imread('ub8.png')
        fig, ax = plt.subplots(figsize=(8,4),dpi=320)
        # fig.figure(figsize=(8, 4))
        # figure(num=None, figsize=(8, 4), dpi=80, facecolor='w', edgecolor='k')
        ax.imshow(im, aspect='auto', zorder=-1, extent=[2.1,4.1,-1,2.1])
        # ax.imshow(im, aspect='auto', extent=(0.4, 0.6, .5, .7), zorder=-1)
        lc = LineCollection(segments, cmap=plt.get_cmap('viridis'))
        lc.set_array(x)
        lc.set_linewidth(3)
        ax.add_collection(lc)
        ax.autoscale_view()
        plt.xlim([2.2, 4])
        plt.xlabel('r ($nm$)', fontsize=18)
        plt.ylabel('g(r)', fontsize=18)
        plt.xticks([2.2, 2.6, 3, 3.4, 3.8],fontsize=12)
        plt.yticks([0, 0.4, 0.8, 1.2, 1.6],fontsize=12)
        plt.ylim([0, 1.8])
        plt.yticks(fontsize=12)
        plt.savefig('../TOC.pdf',  bbox_inches='tight')        
        # plt.show()


