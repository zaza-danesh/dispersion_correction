#!/usr/bin/env python3

#running command: python soup_analysis.py -gmx gmx_seq -make_ndx do -msd do -rotacf do -sasa do
import argparse, os, subprocess
import numpy as np
# import matplotlib
# matplotlib.use('GTKAgg')  #I had to use GTKAgg for this to work, GTK threw errors
# import matplotlib.pyplot as plt  #...  and now do whatever you need...
# from decimal import Decimal


def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-gmx", "--gromacs_command", help="",   type=str,  default="gmx")
    parser.add_argument("-traj_treat_rot_trans", "--traj_treat_rot_trans", help="",   type=str,  default=None)
    parser.add_argument("-traj_treat_protein", "--traj_treat_protein", help="",   type=str,  default=None)
    parser.add_argument("-disre_traj", "--disre_traj", help="",   type=str,  default=None)
    parser.add_argument("-traj_treat_whole_center", "--traj_treat_whole_center", help="",   type=str,  default=None)
    parser.add_argument("-traj_analysis", "--traj_analysis", help="",   type=str,  default=None)
    parser.add_argument("-make_ndx", "--make_ndx", help="",   type=str,  default=None)
    parser.add_argument("-make_ndx_NH", "--make_ndx_NH", help="",   type=str,  default=None)
    parser.add_argument("-rotacf", "--rotacf", help="",   type=str,  default=None)
    parser.add_argument("-S2", "--S2", help="",   type=str,  default=None)
    parser.add_argument("-S2_plot", "--S2_plot", help="",   type=str,  default=None)
    parser.add_argument("-msd", "--msd", help="",   type=str,  default=None)
    parser.add_argument("-msd_water", "--msd_water", help="",   type=str,  default=None)
    parser.add_argument("-repeat", "--repeat", help="",   type=str,  default=None)
    parser.add_argument("-msd_plot", "--msd_plot", help="",   type=str,  default=None)
    parser.add_argument("-gyrate", "--gyrate", help="",   type=str,  default=None)
    parser.add_argument("-gyrate_plot", "--gyrate_plot", help="",   type=str,  default=None)
    parser.add_argument("-hbond_pw", "--hbond_pw", help="",   type=str,  default=None)
    parser.add_argument("-hbond_pp", "--hbond_pp", help="",   type=str,  default=None)
    parser.add_argument("-hbond_pp_plot", "--hbond_pp_plot", help="",   type=str,  default=None)
    parser.add_argument("-hbond_pw_plot", "--hbond_pw_plot", help="",   type=str,  default=None)
    parser.add_argument("-hbond_pw_ac", "--hbond_pw_ac", help="",   type=str,  default=None)
    parser.add_argument("-hbond_pp_ac", "--hbond_pp_ac", help="",   type=str,  default=None)
    parser.add_argument("-hbond_ac_plot", "--hbond_ac_plot", help="",   type=str,  default=None)
    parser.add_argument("-energetics", "--energetics", help="",   type=str,  default=None)
    parser.add_argument("-energetics_plot", "--energetics_plot", help="",   type=str,  default=None)
    parser.add_argument("-mindist", "--mindist", help="",   type=str,  default=None)
    parser.add_argument("-mindist_plot", "--mindist_plot", help="",   type=str,  default=None)
    parser.add_argument("-asso_rate_plot", "--asso_rate_plot", help="",   type=str,  default=None)
    parser.add_argument("-rms", "--rms", help="",   type=str,  default=None)
    parser.add_argument("-rmsf", "--rmsf", help="",   type=str,  default=None)
    parser.add_argument("-sasa", "--sasa", help="",   type=str,  default=None)
    parser.add_argument("-sasa_plot", "--sasa_plot", help="",   type=str,  default=None)
    parser.add_argument("-disre", "--disre", help="",   type=str,  default=None)
    parser.add_argument("-num_viol_plot", "--num_viol_plot", help="",   type=str,  default=None)
    parser.add_argument("-simulation_time", "--simulation_time", help="",   type=str,  default=None)


    args = parser.parse_args()
    return args




if __name__ == '__main__':
    args  = parseArguments()
    traj_directory, md_name = [], []
    gromacs_command = args.gromacs_command
    B=100000
    E=300000

    with open('traj_list.txt') as fn:
        for line in fn:
            line_split = line.split("/")
            traj_directory.append(line.strip())

    ff_labels = []
    for i in range(0, len(traj_directory)):
        ff_labels.append(str(traj_directory[i]))

    if args.make_ndx:
        for i in range(0, len(traj_directory)):
            for j in 1, 2, 4, 8:
                md_name ="ub"+str(j)
                with open("make_index.txt", "w") as fn:
                    for n in range(0, j):
                        fn.write('a '+ str(n*1231+1)+ '-'+ str((n+1)*1231)+"\n")
                        fn.write('name '+ str(19+n)+" chain"+str(n+1)+"\n")
                    fn.write('q\n')
                os.system(gromacs_command+' make_ndx'
                 +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'..tpr'
                 +' -o '+traj_directory[i]+md_name+'_prod.'+str(rep)+'..ndx < make_index.txt >/dev/null 2>/dev/null')
    
    if args.make_ndx_NH:
        with open("make_index_NH.txt", 'w') as fn:
                    fn.write("del 0-24\n")
                    fn.write("a N | a H | a H1 & ! r PRO\n")
                    fn.write("name 0 NH\n")
                    fn.write("q\n")        
        for i in range(0, len(traj_directory)):
            for j in 1, 2, 4, 8:
                md_name = "ub"+str(j)
                os.system(gromacs_command+' make_ndx'
                 +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                 +' -o '+traj_directory[i]+md_name+'_prod_NH.ndx < make_index_NH.txt >/dev/null 2>/dev/null')


                
    # if args.traj_treat:
        # os.system('export GMX_MAXBACKUP=-1')
    if args.traj_treat_whole_center:
        print("\n")
        for i in range(0, len(traj_directory)):
            for j in 1, 2, 4, 8:
                md_name ="ub"+str(j)
                os.system('echo chain1 System | '+args.gromacs_command+' trjconv '
                 +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                 +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                 +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                 +' -pbc mol -ur compact -center'
                 +' -o '+traj_directory[i]+md_name+'_prod_center'+str(rep)+'.gro'
                 +' -b 0 -e 0 >/dev/null 2>/dev/null')
                os.system('echo chain1 System | '+args.gromacs_command+' trjconv '
                 +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                 +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                 +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                 +' -pbc mol -ur compact -center'
                 +' -o '+traj_directory[i]+md_name+'_prod_center'+str(rep)+'.xtc'
                 +' -b '+str(B)+' -e '+str(E))
                 # +' >/dev/null 2>/dev/null')
    
    if args.traj_treat_rot_trans:
        for i in range(0, len(traj_directory)):
            # for j in 1, 2, 4, 8:
            j=1
            md_name ="ub"+str(j)
            os.system('echo Protein Protein | '+args.gromacs_command+' trjconv '
             +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
             +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
             +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
             +' -fit rot+trans'
             +' -o '+traj_directory[i]+md_name+'_prod_rot_trans.xtc'
             +' -b '+str(B)+' -e '+str(E)#)
             +' >/dev/null 2>/dev/null')

        print("traj_treat_rot_trans done")
    
    if args.disre_traj:
        for i in range(0, len(traj_directory)):
            # for j in 1, 2, 4, 8:
            j=1
            md_name ="ub"+str(j)
            os.system('echo System System | '+args.gromacs_command+' trjconv '
             +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
             +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
             +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
             +' -o '+traj_directory[i]+md_name+'_prod_protein.xtc'
             +' -b '+str(B)+' -e '+str(E)#)
             +' >/dev/null 2>/dev/null')

            os.system('echo System System | '+args.gromacs_command+' trjconv '
             +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
             +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
             +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
             +' -o '+traj_directory[i]+md_name+'_prod_protein.gro'
             +' -b '+str(B)+' -e '+str(B)#)
             +' >/dev/null 2>/dev/null')

    if args.traj_treat_protein:
        for i in range(0, len(traj_directory)):
            # for j in 1, 2, 4, 8:
            j=1
            md_name ="ub"+str(j)
            os.system('cp '+traj_directory[i]+md_name+'.top '+traj_directory[i]+md_name+'_disre.top')
            # os.system('echo Protein Protein | '+args.gromacs_command+' trjconv '
            #  +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
            #  +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
            #  +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
            #  +' -o '+traj_directory[i]+md_name+'_prod_protein.xtc'
            #  +' -b '+str(B)+' -e '+str(E)#)
            #  +' >/dev/null 2>/dev/null')


    if args.rotacf:
        for i in range(0, len(traj_directory)):
            # for j in 2, 4, 8:
            j=1
            md_name ="ub"+str(j)
            os.system(args.gromacs_command+' rotacf '
             +' -f '+traj_directory[i]+md_name+'_prod_rot_trans.xtc'
             +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
             +' -n '+traj_directory[i]+md_name+'_prod_NH.ndx'
             +' -noaver -d -P 2'
             +' -o '+traj_directory[i]+md_name+'_rotacf.xvg'
             +' -b '+str(B)+' -e '+str(E)#)
             +' >/dev/null 2>/dev/null')
        print("rotacf done")


    if args.S2:
        for i in range(0, len(traj_directory)):
            # for j in 1, 2, 4, 8:
            j=1
            md_name ="ub"+str(j)
            os.system("./process_rotacf.pl "
             +traj_directory[i]+md_name+'_rotacf.xvg '
             +traj_directory[i]+md_name+'_S2.xvg ')
            if j == 10:
                for c in 1, 19, 23, 24, 25, 31, 37, 38, 53, 72, 74, 76, 76:    
                    c = c-1
                    f = open(traj_directory[i]+md_name+'_S2.xvg', "r")
                    contents = f.readlines()
                    # contents = lines.split()
                    f.close()
                    contents.insert(c, " "+str(c+1)+"   NaN\n")
                    f = open(traj_directory[i]+md_name+'_S2.xvg', "w")
                    contents = "".join(contents)
                    f.write(contents)
                    f.close()
        print("S2 done")

    if args.S2_plot:
        res = list(range(1, 74))
        for i in range(0, len(traj_directory)):
            S2 = []
            j=1
            md_name ="ub"+str(j)                
            with open(traj_directory[i]+md_name+'_S2.xvg') as fn:
                for line in fn:
                    data = line.split()
                    S2.append(float(data[1]))
            plt.plot(res, S2, label=traj_directory[i])

        S2 = []
        with open('S2_exp.xvg') as fn:
                for line in fn:
                    data = line.split()
                    # print line.split()
                    S2.append(float(data[1]))
        plt.plot(res, S2, 'o', label="exp")
        plt.legend(loc=4)
        plt.savefig('S2.pdf', bbox_inches='tight')
        plt.close() 
        print("S2_plot done")
        



    if args.msd:
        L=10 #each 10 ns
        # T=t/L 
        T=200/L  #500 ns 
        # I=(100-percentage)*T/100
        I=100/L
        for i in range(0, len(traj_directory)):
            print(traj_directory[i])
            for j in 1, 2, 4, 8:
                md_name ="ub"+str(j)
                for t in range(I,T):
                    os.system('echo Protein | '+gromacs_command+' msd '
                     +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                     +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                     +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                     +' -o '+traj_directory[i]+md_name+'_msd_'+str(L*(t+1))+'.xvg'
                     + ' -tu ns '+' -ngroup 1'
                     + ' -b '+str(L*t)+' -e '+str(L*(t+1)))
    if args.msd_water:
        for i in range(0, len(traj_directory)):
            for j in 1, 2, 4, 8:
                md_name ="ub"+str(j)
                os.system('echo Water | '+gromacs_command+' msd '
                 +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                 +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                 +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                 +' -o '+traj_directory[i]+md_name+'_msd_water.xvg'                     
                 +' -b '+str(B)+' -e '+str(E))
                     

    if args.msd_plot:
        diff_str = []
        ub=[1, 2 ,4, 8]
        Diff = [[0 for x in range(4)] for y in range(len(traj_directory))] 
        std = [[0 for x in range(4)] for y in range(len(traj_directory))] 
        for i in range(0, len(traj_directory)):
            # print traj_directory[i]
            c=-1
            for j in 1, 2, 4, 8:
                c+=1
                # os.system('less '+traj_directory[i]+'ub'+str(j)+'*msd_100*'+'  | grep "\[" | awk \'{print $5}\' ')
                diff_str.append(subprocess.Popen(['less '+traj_directory[i]+'ub'+str(j)+'_msd_*0.xvg'+'  | grep "\[" | awk \'{print $5}\' '], stdout=subprocess.PIPE, shell=True).communicate()[0].split())
                # print (len(diff_str[-1]))
                diff_flt = [float(ii) for ii in diff_str[-1]]
                # print 'ub', str(j), ' ' , np.mean(diff_flt)
                Diff[i][c] = np.mean(diff_flt)                
                std[i][c] = np.std(diff_flt)                
        for i in range(0, len(traj_directory)):
            # plt.plot(ub,  Diff[i][:]/Diff[i][0], label=traj_directory[i])
            plt.errorbar(ub,  Diff[i][:], std[i][:], label=traj_directory[i])
        
        plt.legend()
        plt.savefig('diff_protein.pdf', bbox_inches='tight')


    if args.gyrate:
        for i in range(0, len(traj_directory)):
            for j in 2, 4, 8:
            # j=1
                md_name ="ub"+str(j)
                for n in range(0, j):
                    os.system('echo chain'+str(n+1)+' | '+gromacs_command+' gyrate'
                     +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                     +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                     +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                     +' -o '+traj_directory[i]+md_name+'_gyrate_'+str(j)+'.'+str(n+1)+'.xvg'
                     +' -b '+str(B)+' -e '+str(E))
                     # +' >/dev/null 2>/dev/null')
    
    if args.gyrate_plot:
        x = list(range(1, len(traj_directory)+1))
        for j in 1, 2, 4, 8:            
            all_y = []
            std = []
            for i in range(0, len(traj_directory)):
                md_name ="ub"+str(j)                
                y=[]
                for n in range(0, j):
                    with open(traj_directory[i]+md_name+'_gyrate_'+str(j)+'.'+str(n+1)+'.xvg') as fn:
                        for line in fn:
                            if line.find("@") < 0 and line.find("#") < 0:
                                data = line.split()
                                y.append(float(data[1]))
                        # print np.mean(g)                    
                all_y.append(np.mean(y))
                std.append(np.std(y)/np.sqrt(j))
            plt.errorbar(x, all_y, std, fmt='o', label=j)
            plt.legend()
            plt.xlim([0, len(traj_directory)+1])
            plt.xticks(x, ff_labels, rotation=75)
            plt.ylabel("gyration radius (nm)")
            plt.savefig('gyration_radius.pdf', bbox_inches='tight')
        plt.close()



    if args.hbond_pp:
        for i in range(0, len(traj_directory)):
            for j in 1, 2, 4, 8:
            # j=2
                md_name ="ub"+str(j)
                os.system('echo Protein Protein | '+gromacs_command+' hbond '
                 +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                 +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                 +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                 +' -num '+traj_directory[i]+md_name+'_hbond_pp.xvg'
                 +' -b '+str(B)+' -e '+str(E))
                 # +' >/dev/null 2>/dev/null')
                

    
    if args.hbond_pp_plot:
        x = list(range(1, len(traj_directory)+1))
        # for j in 1, 2, 4, 8:            
        j=2
        all_y = []
        std = []
        for i in range(0, len(traj_directory)):
            md_name ="ub"+str(j)                
            y=[]
            with open(traj_directory[i]+md_name+'_hbond_pp.xvg') as fn:
                for line in fn:
                    if line.find("@") < 0 and line.find("#") < 0:
                        data = line.split()
                        y.append(float(data[1]))
            all_y.append(np.mean(y))
            std.append(np.std(y))
        # plt.plot(x, all_g, label=j)
        plt.errorbar(x, all_y, std,fmt='o', label=j)
        plt.legend()
        plt.xlim([0, len(traj_directory)+1])
        plt.xticks(x, ff_labels, rotation=75)
        plt.ylabel("HB (protein-protein) (#)")
        plt.savefig('hbond_pp.pdf', bbox_inches='tight')
        plt.close()

    if args.hbond_pw:
        for i in range(0, len(traj_directory)):
            for j in 1, 2, 4, 8:
            # j=2
                md_name ="ub"+str(j)
                os.system('echo Protein Water | '+gromacs_command+' hbond '
                 +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                 +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                 +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                 +' -num '+traj_directory[i]+md_name+'_hbond_pw.xvg'
                 +' -b '+str(B)+' -e '+str(E))
                 # +' >/dev/null 2>/dev/null')
                

    if args.hbond_pw_plot:
        x = list(range(1, len(traj_directory)+1))
        # for j in 1, 2, 4, 8:            
        j=2
        all_y = []
        std = []
        for i in range(0, len(traj_directory)):
            md_name ="ub"+str(j)                
            y=[]
            with open(traj_directory[i]+md_name+'_hbond_pw.xvg') as fn:
                for line in fn:
                    if line.find("@") < 0 and line.find("#") < 0:
                        data = line.split()
                        y.append(float(data[1]))
            all_y.append(np.mean(y))
            std.append(np.std(y))
        # plt.plot(x, all_g, label=j)
        plt.errorbar(x, all_y, std, fmt='o', label=j)
        plt.legend()
        plt.xlim([0, len(traj_directory)+1])
        plt.xticks(x, ff_labels, rotation=75)
        plt.ylabel("HB (protein-water) (#)")
        plt.savefig('hbond_pw.pdf', bbox_inches='tight')



    if args.hbond_pw_ac:
        # for j in 1, 2, 4, 8:
        j=2
        for i in range(0, len(traj_directory)):
            md_name ="ub"+str(j)
            os.system('echo chain1 Water | '+gromacs_command+' hbond '
             +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
             +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
             +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
             +' -ac '+traj_directory[i]+md_name+'_hbond_pw_ac.xvg'
             +' -b '+str(B)+' -e '+str(E)#)
             +' >> pw_ac.txt')
                


    if args.hbond_pp_ac:
        # for j in 1, 2, 4, 8:
        j=2
        for i in range(0, len(traj_directory)):
            md_name ="ub"+str(j)
            os.system('echo chain1 chain2 | '+gromacs_command+' hbond '
             +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
             +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
             +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
             +' -ac '+traj_directory[i]+md_name+'_hbond_pp_ac.xvg'
             +' -b '+str(B)+' -e '+str(E)#)
             +' >> pp_ac.txt')


    if args.hbond_ac_plot:
        ub=[1, 2 ,4, 8]
        pp_forward = subprocess.Popen(['cat pp_ac.txt | grep Forward | awk \'{print $3}\' '],stdout=subprocess.PIPE, shell=True).communicate()[0].split()
        pp_ac = np.split(np.asarray(pp_forward), len(traj_directory))
        
        pw_forward = subprocess.Popen(['cat pw_ac.txt | grep Forward | awk \'{print $3}\''],stdout=subprocess.PIPE, shell=True).communicate()[0].split()
        pw_ac = np.split(np.asarray(pw_forward), len(traj_directory))
        
        for i in range(0, len(traj_directory)):
           plt.plot(ub,  pp_ac[i], fmt='o', label=traj_directory[i])

        plt.legend(loc=1)
        plt.savefig('pp_ac.pdf', bbox_inches='tight')

        # for i in range(0, len(traj_directory)):
        #    plt.plot(ub,  pw_ac[i], 'o--', label=traj_directory[i])

        # plt.legend(loc=1)
        # plt.savefig('pw_ac.pdf', bbox_inches='tight')

    if args.energetics:
        for i in range(0, len(traj_directory)):
            for j in 1, 2, 4, 8:
                md_name ="ub"+str(j)
                os.system('echo Density | '+gromacs_command+' energy '
                 +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.edr'
                 +' -o '+traj_directory[i]+md_name+'_density.xvg'
                 +' -b '+str(B)+' -e '+str(E))
                 # +' >/dev/null 2>/dev/null')
                os.system('echo Potential | '+gromacs_command+' energy '
                 +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.edr'
                 +' -o '+traj_directory[i]+md_name+'_potential.xvg'
                 +' -b '+str(B)+' -e '+str(E))
                 # +' >/dev/null 2>/dev/null')
                os.system('echo Total-Energy | '+gromacs_command+' energy '
                 +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.edr'
                 +' -o '+traj_directory[i]+md_name+'_total_energy.xvg'
                 +' -b '+str(B)+' -e '+str(E))
                 # +' >/dev/null 2>/dev/null')
                 

    if args.energetics_plot:        
        x = list(range(1, len(traj_directory)+1))
        for j in 1, 2, 4, 8:            
            md_name ="ub"+str(j)                
            all_y = []
            std = []
            for i in range(0, len(traj_directory)):
                y=[]
                with open(traj_directory[i]+md_name+'_density.xvg') as fn:
                    for line in fn:
                        if line.find("@") < 0 and line.find("#") < 0:
                            data = line.split()
                            y.append(float(data[1]))
                            # print y[-1]
                all_y.append(np.mean(y))
                # print all_y[-1]
                std.append(np.std(y))
            # plt.plot(x, all_g, label=j)
            # print all_y
            if j==1:
                plt.subplot2grid((8,1), (0,0), rowspan=2)
                plt.ylim([1002, 1012])
                plt.errorbar(x, all_y, std,fmt='o', label=j)
                plt.tick_params(axis='y', which='major', labelsize=10)
                plt.xticks(x, [])
                plt.legend(fontsize=10)
            elif j==2:
                plt.subplot2grid((8,1), (2,0), rowspan=2)
                plt.errorbar(x, all_y, std,fmt='o', label=j, color='r')
                plt.ylim([1009, 1019])
                plt.tick_params(axis='y', which='major', labelsize=10)
                plt.ylabel("(kg/m^3)", fontsize=16)
                plt.legend(prop={'size': 10})
                plt.xticks(x, [])
            elif j==4:
                plt.subplot2grid((8,1), (4,0),rowspan=2)                
                plt.ylim([1026, 1036])
                plt.errorbar(x, all_y, std,fmt='o', label=j, color='g')
                plt.ylabel("density", fontsize=16)
                plt.tick_params(axis='y', which='major', labelsize=10)
                plt.legend(fontsize=10)
                plt.xticks(x, [])
            else:
                plt.subplot2grid((8,1), (6,0), rowspan=2)
                plt.ylim([1055, 1068])
                plt.errorbar(x, all_y, std,fmt='o', label=j, color='m')
                # plt.xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5], ff_labels, rotation=45)
                plt.xticks(x, ff_labels, rotation=75)
                plt.tick_params(axis='y', which='major', labelsize=10)
                plt.legend(fontsize=10)
                # plt.xticks(x, [])

            # plt.legend()
            plt.xlim([0, len(traj_directory)+2])

        plt.savefig('density.pdf', bbox_inches='tight')
        plt.close() 



        for j in 1, 2, 4, 8:            
            md_name ="ub"+str(j)                
            all_y = []
            std = []
            for i in range(0, len(traj_directory)):
                y=[]
                with open(traj_directory[i]+md_name+'_potential.xvg') as fn:
                    for line in fn:
                        if line.find("@") < 0 and line.find("#") < 0:
                            data = line.split()
                            y.append(float(data[1]))
                            # print y[-1]
                all_y.append(np.mean(y))
                # print all_y[-1]
                std.append(np.std(y))
            # plt.plot(x, all_g, label=j)
            # print all_y
            if j==1:
                plt.subplot2grid((8,1), (0,0), rowspan=2)
                # plt.ylim([1000, 1015])
                plt.errorbar(x, all_y, std,fmt='o', label=j)
                plt.tick_params(axis='y', which='major', labelsize=10)
                plt.xticks(x, [])
                plt.legend(fontsize=10)
            elif j==2:
                plt.subplot2grid((8,1), (2,0), rowspan=2)
                plt.errorbar(x, all_y, std,fmt='o', label=j, color='r')
                # plt.ylim([1005, 1020])
                plt.tick_params(axis='y', which='major', labelsize=10)
                plt.ylabel("(kJ/mol)", fontsize=16)
                plt.legend(prop={'size': 10})
                plt.xticks(x, [])
            elif j==4:
                plt.subplot2grid((8,1), (4,0),rowspan=2)                
                # plt.ylim([1025, 1040])
                plt.errorbar(x, all_y, std,fmt='o', label=j, color='g')
                plt.ylabel("potential", fontsize=16)
                plt.tick_params(axis='y', which='major', labelsize=10)
                plt.legend(fontsize=10)
                plt.xticks(x, [])
            else:
                plt.subplot2grid((8,1), (6,0), rowspan=2)
                # plt.ylim([1055, 1070])
                plt.errorbar(x, all_y, std,fmt='o', label=j, color='m')
                # plt.xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5], ff_labels, rotation=45)
                plt.xticks(x, ff_labels, rotation=75)
                plt.tick_params(axis='y', which='major', labelsize=10)
                plt.legend(fontsize=10)
                # plt.xticks(x, [])

            # plt.legend()
            plt.xlim([0, len(traj_directory)+2])

        plt.savefig('potential.pdf', bbox_inches='tight')
        plt.close()


    if args.mindist:
        for i in range(0, len(traj_directory)):
            # for j in 1, 2, 4, 8:
            j=2
            md_name ="ub"+str(j)
            os.system('echo chain1 chain2 | '+gromacs_command+' mindist '
                +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                +' -od '+traj_directory[i]+md_name+'_mindist.xvg'
                +' -on '+traj_directory[i]+md_name+'_numcont.xvg'
                +' -o '+traj_directory[i]+md_name+'_atom-pair.out'
                # +' -ox '+traj_directory[i]+md_name+'_mindist.xtc'
                +' -or '+traj_directory[i]+md_name+'_mindistres.xvg'
                +'  -b '+str(B)+' -e '+str(E))
                # +' >/dev/null 2>/dev/null')



    if args.mindist_plot:
        x = list(range(1, len(traj_directory)+1))
        j=2
        all_y = []
        std = []
        for i in range(0, len(traj_directory)):
            md_name ="ub"+str(j)                
            y=[]
            with open(traj_directory[i]+md_name+'_mindist.xvg') as fn:
                for line in fn:
                    if line.find("@") < 0 and line.find("#") < 0:
                        data = line.split()
                        y.append(float(data[1]))
            all_y.append(np.mean(y))
            std.append(np.std(y))
        # plt.plot(x, all_g, label=j)
        plt.errorbar(x, all_y, std, fmt='o', label=j)
        plt.legend()
        plt.xlim([0, len(traj_directory)+1])
        plt.xticks(x, ff_labels, rotation=75)
        plt.ylabel("min distance (nm)")
        plt.savefig('mindist.pdf', bbox_inches='tight')
        plt.close()

        all_y = []
        std = []
        for i in range(0, len(traj_directory)):
            md_name ="ub"+str(j)                
            y=[]
            with open(traj_directory[i]+md_name+'_numcont.xvg') as fn:
                for line in fn:
                    if line.find("@") < 0 and line.find("#") < 0:
                        data = line.split()
                        y.append(float(data[1]))
            all_y.append(np.mean(y))
            std.append(np.std(y))
        # plt.plot(x, all_g, label=j)
        plt.errorbar(x, all_y, std, fmt='o', label=j)
        plt.legend()
        plt.xlim([0, len(traj_directory)+1])
        plt.xticks(x, ff_labels, rotation=75)
        plt.ylabel("number contact (dist < 0.6 nm )")
        plt.savefig('numcont.pdf', bbox_inches='tight')


    if args.asso_rate_plot:
        MINDIST = 0.25
        x = list(range(1, len(traj_directory)+1))
        j=2
        all_y = []
        ste = []
        for i in range(0, len(traj_directory)):
            md_name ="ub"+str(j)                
            y=[0]
            # T=[0]
            t_temp=0
            with open(traj_directory[i]+md_name+'_mindist.xvg') as fn:
                for line in fn:
                    if line.find("@") < 0 and line.find("#") < 0:
                        data = line.split()
                        t = t_temp
                        if float(data[1]) < MINDIST:
                            t_temp = float(data[0])
                            # print float(data[1]), t 
                        else:
                            if y[-1] != t:
                                y.append(t)
                                # print T[-1]

            asso_time = np.diff(y)
            print asso_time
            all_y.append(np.mean(asso_time))
            ste.append(np.std(asso_time)/np.sqrt(len(asso_time)))
            print all_y[-1], ste[-1]

        plt.errorbar(x, all_y, ste, fmt='o', label=j)
        plt.legend()
        plt.xlim([0, len(traj_directory)+1])
        plt.xticks(x, ff_labels, rotation='vertical')
        plt.ylabel("average contact time")
        plt.savefig('asso_time.pdf', bbox_inches='tight')
        plt.close()

    if args.rms:
        for i in range(0, len(traj_directory)):
            for j in 1, 2, 4, 8:
                md_name ="ub"+str(j)
                os.system('echo chain1 chain1 | '+gromacs_command+' rms '
                    +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                    +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                    +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                    +' -o '+traj_directory[i]+md_name+'_rms.xvg'
                    +' -b '+str(B)+' -e '+str(E))
                    # +' >/dev/null 2>/dev/null')

    if args.rmsf:
        for i in range(0, len(traj_directory)):
            for j in 1, 2, 4, 8:
                md_name ="ub"+str(j)
                os.system('echo chain1 chain1 | '+gromacs_command+' rmsf '
                    +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                    +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                    +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                    +' -o '+traj_directory[i]+md_name+'_rmsf.xvg'
                    +' -b '+str(B)+' -e '+str(E))
                    # +' >/dev/null 2>/dev/null')

    if args.sasa:
        for i in range(0, len(traj_directory)):
            for j in 1, 2, 4, 8:
                md_name ="ub"+str(j)
                os.system('echo Protein | '+gromacs_command+' sasa '
                    +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                    +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                    +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                    +' -o '+traj_directory[i]+md_name+'_sasa.xvg'
                    +' -or '+traj_directory[i]+md_name+'_sasa_res.xvg'
                    +' -tv '+traj_directory[i]+md_name+'_sasa_tv.xvg'
                    +' -b '+str(B)+' -e '+str(E))
                    # +' >/dev/null 2>/dev/null')

    if args.sasa_plot:
        x = list(range(1, len(traj_directory)+1))
        for j in 1, 2, 4, 8:            
            all_y = []
            std = []
            for i in range(0, len(traj_directory)):
                md_name ="ub"+str(j)                
                y=[]
                for n in range(0, j):
                    with open(traj_directory[i]+md_name+'_sasa.xvg') as fn:
                        for line in fn:
                            if line.find("@") < 0 and line.find("#") < 0:
                                data = line.split()
                                y.append(float(data[1]))
                        # print np.mean(g)                    
                all_y.append(np.mean(y))
                std.append(np.std(y)/np.sqrt(j))
            plt.errorbar(x, all_y, std, fmt='o', label=j)
            plt.legend()
            plt.xlim([0, len(traj_directory)+1])
            plt.xticks(x, ff_labels, rotation='vertical')
            plt.ylabel("sasa (nm^2)")
            plt.savefig('sasa.pdf', bbox_inches='tight')
        plt.close()



    if args.disre:
        for i in range(0, len(traj_directory)):
            # for j in 1, 2, 4, 8:
            j=1
            md_name ="ub"+str(j)
            os.system(args.gromacs_command+' disre '
                +' -f '+traj_directory[i]+md_name+'_prod_protein.xtc'
                +' -s '+traj_directory[i]+md_name+'_disre.tpr'
                +' -ds '+traj_directory[i]+md_name+'_drsum.xvg'
                +' -da '+traj_directory[i]+md_name+'_draver.xvg'
                +' -dn '+traj_directory[i]+md_name+'_drnum.xvg'
                +' -dm '+traj_directory[i]+md_name+'_drmax.xvg'
                +' -dr '+traj_directory[i]+md_name+'_restr.xvg'
                +' -l '+traj_directory[i]+md_name+'_disres.log'
                +' -q '+traj_directory[i]+md_name+'_viol.pdb')


    if args.num_viol_plot:
        all_y = []
        x = list(range(1, len(traj_directory)+1))        
        for i in range(0, len(traj_directory)):
            temp = []
            # for j in 1, 2, 4, 8:
            j=1
            md_name ="ub"+str(j)
            temp = (subprocess.Popen(['cat '+traj_directory[i]+md_name+'_disres.log '+' | grep "Number" | awk \'{print $5}\' '], stdout=subprocess.PIPE, shell=True).communicate()[0].split())
            
            all_y.append(float(temp[2].split("/")[0]))
        total = float(temp[2].split("/")[1])
        all_y[:] = [el / total for el in all_y]
        
        # print all_y
        print len(all_y)
        print len(x)
        plt.plot(x, all_y, 'o', label=j)
        # plt.legend()
        plt.xlim([0, len(traj_directory)+1])
        plt.xticks(x, ff_labels, rotation=75)
        plt.ylabel("number of violated restraints (%)")
        plt.savefig('num_viol.pdf', bbox_inches='tight')
        plt.close()

        # print float(all_y)/float(total)

    if args.simulation_time:
        for i in range(0, len(traj_directory)):
            print()
            print(traj_directory[i])
            for j in 1, 2, 4, 8:
                for rep in 1, 2, 3:
                    t = subprocess.Popen(['tail -n 400 '+traj_directory[i]+'ub'+str(j)+'_prod.'+str(rep)+'.log | grep "DD  step" | tail -n 1| awk \'{print $3}\' '], stdout=subprocess.PIPE, shell=True).communicate()[0].split()
                    print("ub",j,".",rep,":",float(t[0])*0.002)
                    # print 'ub'+str(j)+'_prod.'+str(rep)+float(t[0])*0.002
                # print float(t)*0.002
