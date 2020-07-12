#!/usr/bin/env python3
import argparse, os, subprocess

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-gmx", "--gromacs_command", help="",   type=str,  default="gmx")
    parser.add_argument("-traj_list", "--traj_list", help="",   type=str,  default=None)
    parser.add_argument("-make_ndx", "--make_ndx", help="",   type=str,  default=None)
    parser.add_argument("-energetics", "--energetics", help="",   type=str,  default=None)
    parser.add_argument("-mindist", "--mindist", help="",   type=str,  default=None)
    parser.add_argument("-sasa", "--sasa", help="",   type=str,  default=None)
    parser.add_argument("-gyrate", "--gyrate", help="",   type=str,  default=None)
    parser.add_argument("-rms", "--rms", help="",   type=str,  default=None)
    parser.add_argument("-rmsdist", "--rmsdist", help="",   type=str,  default=None)
    parser.add_argument("-rmsf", "--rmsf", help="",   type=str,  default=None)
    parser.add_argument("-treat_traj", "--treat_traj", help="",   type=str,  default=None)
    parser.add_argument("-hbond_pp", "--hbond_pp", help="",   type=str,  default=None)
    parser.add_argument("-rdf", "--rdf", help="",   type=str,  default=None)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args  = parseArguments()
    traj_directory, md_name = [], []
    gromacs_command = args.gromacs_command
    B=100000
    E=278000

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
                for rep in 1, 2, 3:
                    with open("make_index.txt", "w") as fn:
                        for n in range(0, j):
                            fn.write('a '+ str(n*1231+1)+ '-'+ str((n+1)*1231)+"\n")
                            fn.write('name '+ str(19+n)+" chain"+str(n+1)+"\n")
                        fn.write('q\n')
                    os.system(gromacs_command+' make_ndx'
                     +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                     +' -o '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx < make_index.txt >/dev/null 2>/dev/null')

    if args.energetics:
        for i in range(0, len(traj_directory)):
            for j in 1, 2, 4, 8:
                md_name ="ub"+str(j)
                for rep in 1, 2, 3:
                    # os.system('echo Density | '+gromacs_command+' energy '
                    #  +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.edr'
                    #  +' -o '+traj_directory[i]+md_name+'_density.'+str(rep)+'.xvg'
                    #  # +' -b '+str(B)+' -e '+str(E)
                    #  +' >/dev/null 2>/dev/null')
                    # os.system('echo Potential | '+gromacs_command+' energy '
                    #  +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.edr'
                    #  +' -o '+traj_directory[i]+md_name+'_potential.'+str(rep)+'.xvg'
                    #  # +' -b '+str(B)+' -e '+str(E)
                    #  +' >/dev/null 2>/dev/null')
                    # os.system('echo Total-Energy | '+gromacs_command+' energy '
                    #  +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.edr'
                    #  +' -o '+traj_directory[i]+md_name+'_total_energy.'+str(rep)+'.xvg'
                    #  # +' -b '+str(B)+' -e '+str(E)
                    #  +' >/dev/null 2>/dev/null')
                    os.system(gromacs_command+' energy '
                     +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.edr'
                     +' -o '+traj_directory[i]+md_name+'_lj.'+str(rep)+'.xvg'
                     +' -b '+str(B)+' -e '+str(E)+' -sum'
                     +' < lj.input >/dev/null 2>/dev/null')

    if args.mindist:
        for i in range(0, len(traj_directory)):
            j=2
            md_name ="ub"+str(j)
            for rep in 1, 2, 3:
                os.system('echo chain1 chain2 | '+gromacs_command+' mindist '
                    +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                    +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                    +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                    +' -od '+traj_directory[i]+md_name+'_mindist.'+str(rep)+'.xvg'
                    # +' -on '+traj_directory[i]+md_name+'_numcont.'+str(rep)+'.xvg'
                    # +' -o '+traj_directory[i]+md_name+'_atom-pair.'+str(rep)+'.out'
                    # +' -ox '+traj_directory[i]+md_name+'_mindist.xtc'
                    # +' -or '+traj_directory[i]+md_name+'_mindistres.'+str(rep)+'.xvg'
                    # +'  -b '+str(B)+' -e '+str(E)
                    +' >/dev/null 2>/dev/null')

    if args.sasa:
        for i in range(0, len(traj_directory)):
            j=2
            for chain in 1, 2:
                md_name ="ub"+str(j)
                for rep in 1, 2, 3:
                    os.system('echo chain'+ str(chain) +' | '+gromacs_command+' sasa '
                        +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                        +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                        +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                        +' -o '+traj_directory[i]+md_name+'.'+str(chain)+'_sasa.'+str(rep)+'.xvg'
                        +' -or '+traj_directory[i]+md_name+'.'+str(chain)+'_sasa_res.'+str(rep)+'.xvg'
                        +' -tv '+traj_directory[i]+md_name+'.'+str(chain)+'_sasa_tv.'+str(rep)+'.xvg'
                        +' -b '+str(B)+' -e '+str(E)
                        +' >/dev/null 2>/dev/null')

    if args.gyrate:
        for i in range(0, len(traj_directory)):
            j=1
            # for j in 2, 4, 8:
            md_name ="ub"+str(j)
            for rep in 1, 2, 3:
                # for chain in range(1, j+1):
                    # print("chain"+str(chain))
                    os.system('echo Protein | '+gromacs_command+' gyrate '
                    # os.system('echo chain'+str(chain)+' | '+gromacs_command+' gyrate '
                        +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                        +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                        +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                        # +' -o '+traj_directory[i]+md_name+'_chain'+str(chain)+'_gyrate.'+str(rep)+'.xvg'
                        +' -o '+traj_directory[i]+md_name+'_gyrate.'+str(rep)+'.xvg'
                        # +' -dt 1000 ')
                        # +' -b '+str(B)+' -e '+str(E) )
                        +' >/dev/null 2>/dev/null')

    if args.rms:
        for i in range(0, len(traj_directory)):
            j=1
            md_name ="ub"+str(j)
            for rep in 1, 2, 3:
                os.system('echo 4 4 | '+gromacs_command+' rms '
                    # +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                    +' -f '+traj_directory[i]+md_name+'_prod_whole_center_rot+trans.'+str(rep)+'.xtc'
                    +' -s '+traj_directory[i]+md_name+'_analyze.'+str(rep)+'.tpr'
                    # +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                    +' -o '+traj_directory[i]+md_name+'_rms.'+str(rep)+'.xvg'
                    # +' -b '+str(B)+' -e '+str(E)
                    +' >/dev/null 2>/dev/null')


    if args.rmsdist:
        for i in range(0, len(traj_directory)):
            j=1
            md_name ="ub"+str(j)
            for rep in 1, 2, 3:
                if args.treat_traj:
                    os.system('echo Protein Protein | '+gromacs_command+' trjconv '
                        +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                        +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                        +' -pbc whole -ur compact -center'
                        +' -dt 1000'
                        +' -o '+traj_directory[i]+md_name+'_prod_whole_center.'+str(rep)+'.xtc')
                    os.system('echo Protein | '+gromacs_command+' trjconv '
                        +' -f '+traj_directory[i]+md_name+'_prod_whole_center.'+str(rep)+'.xtc'
                        +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                        +' -b 100000 -e 100000'
                        +' -o '+traj_directory[i]+md_name+'_prod_whole_center.'+str(rep)+'.gro')
                    # os.system('echo Protein Protein | '+gromacs_command+' trjconv '
                    #     +' -f '+traj_directory[i]+md_name+'_prod_whole_center.'+str(rep)+'.xtc'
                    #     +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                    #     +' -fit rot+trans'
                    #     +' -o '+traj_directory[i]+md_name+'_prod_whole_center_rot+trans.'+str(rep)+'.xtc')                    
                os.system('echo 4 4 | '+gromacs_command+' rmsdist '
                    # +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                    +' -f '+traj_directory[i]+md_name+'_prod_whole_center.'+str(rep)+'.xtc'
                    +' -s '+traj_directory[i]+md_name+'_prod_whole_center.'+str(rep)+'.gro'
                    # +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                    +' -o '+traj_directory[i]+md_name+'_rmsdist_fromgro.'+str(rep)+'.xvg'
                    # +' -dt 100 '
                    +' -b '+str(B)+ ' -e 200000 '
                    +' >/dev/null 2>/dev/null')
    if args.rmsf:
        for i in range(0, len(traj_directory)):
            j=1
            md_name ="ub"+str(j)
            for rep in 1, 2, 3:
                if args.treat_traj:
                    os.system('echo Protein Protein | '+gromacs_command+' trjconv '
                        +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                        +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                        +' -pbc whole -ur compact -center'
                        +' -dt 1000'
                        +' -o '+traj_directory[i]+md_name+'_prod_whole_center.'+str(rep)+'.xtc')
                    os.system('echo Protein Protein | '+gromacs_command+' trjconv '
                        +' -f '+traj_directory[i]+md_name+'_prod_whole_center.'+str(rep)+'.xtc'
                        +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                        +' -fit rot+trans'
                        +' -o '+traj_directory[i]+md_name+'_prod_whole_center_rot+trans.'+str(rep)+'.xtc')

                os.system('echo Protein | '+gromacs_command+' trjconv '
                    +' -f '+traj_directory[i]+md_name+'_prod_whole_center_rot+trans.'+str(rep)+'.xtc'
                    +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                    +' -o '+traj_directory[i]+md_name+'_prod_whole_center_rot+trans.'+str(rep)+'.gro'
                    +' -b '+str(B)+' -e '+str(B))
                os.system('echo Backbone | '+gromacs_command+' rmsf '
                    +' -f '+traj_directory[i]+md_name+'_prod_whole_center_rot+trans.'+str(rep)+'.xtc'
                    +' -s '+traj_directory[i]+md_name+'_prod_whole_center_rot+trans.'+str(rep)+'.gro'
                    +' -o '+traj_directory[i]+md_name+'_rmsf.'+str(rep)+'.xvg'
                    +' -res '
                    +' -b '+str(B))#+' -e '+str(E))
                    # +' >/dev/null 2>/dev/null')


    if args.rdf:
        for i in range(0, len(traj_directory)):
            for j in 2, 4, 8:
                md_name="ub"+str(j)
                for rep in 1, 2, 3:
                    # os.system('echo Protein | '+ gromacs_command+' trjconv '
                    #     +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                    #     +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                    #     +' -o '+traj_directory[i]+md_name+'_prot.'+str(rep)+'.xtc'
                    #     # +' -tu ns -b 100 -e 300 '
                    #     +' -tu ns -b 100  '
                    #     +' >/dev/null 2>/dev/null ')
                    os.system(gromacs_command+' rdf '
                        +' -f '+traj_directory[i]+md_name+'_prot.'+str(rep)+'.xtc'
                        +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                        +' -o '+traj_directory[i]+md_name+'_rdf_b100_e278ns.'+str(rep)+'.xvg'
                        # +' -cn '+traj_directory[i]+md_name+'_rdf-cn.'+str(rep)+'.xvg'
                        +' --selrpos mol_com --seltype whole_mol_com -sel Protein -ref Protein'
                        # +' -tu ns -b 100 -e 300 -dt 0.01 '
                        +' -tu ns -b 100 -e 278 -dt 0.01 '
                        +' >/dev/null 2>/dev/null ')

    if args.hbond_pp:
        for i in range(0, len(traj_directory)):
            j=2
            md_name ="ub"+str(j)
            for rep in 1, 2, 3:
                os.system('echo chain1 chain2 | '+gromacs_command+' hbond '
                    +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                    +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                    +' -n '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.ndx'
                    +' -num '+traj_directory[i]+md_name+'_hb_num_pp.'+str(rep)+'.xvg'
                    +' -b '+str(B)+' -e '+str(E)
                    +' >/dev/null 2>/dev/null')


    if args.trjconv:
        # gmx trjconv -f ffam-ws/ub8_prod.1.xtc -s ffam-ws/ub8_prod.1.tpr -o ub8_prot.1 -tu ns -b 100 -e 300
        for i in range(0, len(traj_directory)):
            for j in 2, 4, 8:
                md_name="ub"+str(j)
                for rep in 1, 2, 3:
                    os.system('echo Protein | '+ gromacs_command+' trjconv '
                        +' -f '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.xtc'
                        +' -s '+traj_directory[i]+md_name+'_prod.'+str(rep)+'.tpr'
                        +' -o '+traj_directory[i]+md_name+'_prot.'+str(rep)+'.xtc'
                        +' -tu ns -b 100 -e 300 '
                        +' >/dev/null 2>/dev/null ')