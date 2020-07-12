import argparse, os


def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-gmx", "--gromacs_command" , help="",   type=str,  default="gmx")
    parser.add_argument("-make_top", "--make_top" , help="",   type=str,  default=None)
    parser.add_argument("-wd", "--wd" , help="working directory",   type=str,  default=None)
    parser.add_argument("-solvate", "--solvate" , help="do solvate",   type=str,  default=None)
    parser.add_argument("-add_salt", "--add_salt" , help="do add_salt",   type=str,  default=None)
    parser.add_argument("-molar", "--molar" , help="ionic strength",   type=str,  default=None)
    parser.add_argument("-CATION", "--CATION" , help="positive ion name",   type=str,  default=None)
    parser.add_argument("-ANION", "--ANION" , help="negative ion name",   type=str,  default=None)
    parser.add_argument("-em", "--em" , help="energy minimization step",   type=str,  default=None)
    parser.add_argument("-nvt", "--nvt" , help="nvt step",   type=str,  default=None)
    parser.add_argument("-npt", "--npt" , help="npt step",   type=str,  default=None)
    parser.add_argument("-prod", "--prod" , help="production step",   type=str,  default=None)
    parser.add_argument("-analyze", "--analyze" , help="production step",   type=str,  default=None)
    parser.add_argument("-rerun_energy", "--rerun_energy" , help="production step",   type=str,  default=None)
    parser.add_argument("-disre", "--disre" , help="distance restraint",   type=str,  default=None)


    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args  = parseArguments()
    mdp_index=""    
    if args.wd == "ffam-cf00":
        forcefield="amber99sb-ildn"
        ch_am = "amber" 
    elif args.wd == "ffam-cf01":
        forcefield="amber99sbildn_pps1"
        ch_am = "amber" 
    elif args.wd == "ffam-cf02":
        forcefield="amber99sbildn_pps2"    
        ch_am = "amber" 
    elif args.wd == "ffam-cf04":
        forcefield="amber99sbildn_pps4"    
        ch_am = "amber" 
    elif args.wd == "ffam-ws":
        forcefield="amber99sbws"
        ch_am = "amber" 
    elif args.wd == "ffam03-ws":
        forcefield="amber03ws"
        ch_am = "amber"
    elif args.wd == "ffam-pme00":
        forcefield="amber99sb-ildn"
        mdp_index="_vdwpme"
        ch_am = "amber" 
    elif args.wd == "ffam-pme01":
        forcefield="amber99sbildn_pps1"
        mdp_index="_vdwpme"
        ch_am = "amber" 
    elif args.wd == "ffam-pme02":
        forcefield="amber99sbildn_pps2"
        mdp_index="_vdwpme"
        ch_am = "amber" 
    elif args.wd == "ffam-pme04":
        forcefield="amber99sbildn_pps4"
        mdp_index="_vdwpme"
        ch_am = "amber" 
    elif args.wd == "ffch-pme01":
        forcefield="charmm27_pps1"
        mdp_index="_vdwpme"
        ch_am = "charmm" 
    elif args.wd == "ffch-pme02":
        forcefield="charmm27_pps2"
        mdp_index="_vdwpme"
        ch_am = "charmm" 
    elif args.wd == "ffch-pme04":
        forcefield="charmm27_pps4"
        mdp_index="_vdwpme"
        ch_am = "charmm" 
    elif args.wd == "ffch-pme00":
        forcefield="charmm27"
        mdp_index="_vdwpme"
        ch_am = "charmm" 
    else:
        print("no such forcefield's parameter prepared")
    # os.system('cp -r '+forcefield+'.ff'+' '+args.wd )
    if ch_am == "amber":
        os.system('cp CONF/'+ch_am+'/* CONF/')
    elif ch_am == "charmm":
        os.system('cp CONF/'+ch_am+'/* CONF/')
    else:
        print("error in forcfield name")

    if args.make_top :
        if args.wd == "ffam-ws" or args.wd == "ffam03-ws":
            for i in 1, 2, 4, 8:
                with open(args.wd+'/ub'+str(i)+'.top', 'w') as fn:
                    fn.write('#include "'+forcefield+'.ff/forcefield.itp" \n')
                    fn.write(' \n')
                    fn.write('#include "../CONF/ubiquitin.itp" \n')
                    fn.write(' \n')
                    fn.write('#include "'+forcefield+'.ff/tip4p2005s.itp" \n')
                    fn.write('#ifdef POSRES_WATER \n')
                    fn.write('[ position_restraints ] \n')
                    fn.write('   1    1       1000       1000       1000 \n')
                    fn.write('#endif \n')
                    fn.write(' \n')
                    fn.write('#include "'+forcefield+'.ff/ions.itp" \n')
                    fn.write(' \n')
                    fn.write(' \n')
                    fn.write('[ system ] \n')
                    fn.write('ub'+str(i)+' \n')
                    fn.write(' \n')
                    fn.write('[ molecules ] \n')
                    fn.write('Ubiquitin\t\t'+str(i)+' \n')

        else:
            for i in 1, 2, 4, 8:
                with open(args.wd+'/ub'+str(i)+'.top', 'w') as fn:
                    fn.write('#include "'+forcefield+'.ff/forcefield.itp" \n')
                    fn.write('#include "../CONF/tip4p2005_types.itp"\n')                
                    fn.write(' \n')
                    fn.write('#include "../CONF/ubiquitin.itp" \n')
                    fn.write(' \n')
                    fn.write('#include "../CONF/tip4p2005.itp"\n')
                    fn.write('#ifdef POSRES_WATER \n')
                    fn.write('[ position_restraints ] \n')
                    fn.write('   1    1       1000       1000       1000 \n')
                    fn.write('#endif \n')
                    fn.write(' \n')
                    fn.write('#include "'+forcefield+'.ff/ions.itp" \n')
                    fn.write(' \n')
                    fn.write(' \n')
                    fn.write('[ system ] \n')
                    fn.write('ub'+str(i)+' \n')
                    fn.write(' \n')
                    fn.write('[ molecules ] \n')
                    fn.write('Ubiquitin\t\t'+str(i)+' \n')
    

    if args.solvate:
        # for i in 1, 2, 4, 8:
        i=1
        os.system('cp CONF/'+ch_am+'/box_ub'+str(i)+'.gro '+args.wd)
        os.system(args.gromacs_command+' solvate -cp '+args.wd+'/box_ub'+str(i)+'.gro -cs CONF/tip4p2005.gro -p '+args.wd+'/ub'+str(i)+'.top -o '+args.wd+'/ub'+str(i)+'_solv.gro ')

    if args.add_salt:
        # for i in 1, 2, 4, 8:
        i=1
        os.system(args.gromacs_command+
            ' grompp -f MDP/em'+mdp_index+'.mdp -c '+args.wd+'/ub'+str(i)+'_solv.gro -p '+args.wd+'/ub'+str(i)+'.top -o '+args.wd+'/ub'+str(i)+'_ions.tpr')
        os.system('echo SOL | '+args.gromacs_command
            +' genion -s '+args.wd+'/ub'+str(i)+'_ions.tpr -p '+args.wd+'/ub'+str(i)+'.top -neutral -conc 0.15 -pname NA -nname CL  -o '+args.wd+'/ub'+str(i)+'_is.gro')
        os.system(args.gromacs_command
            +' grompp -f MDP/em'+mdp_index+'.mdp -c '+args.wd+'/ub'+str(i)+'_is.gro -p '+args.wd+'/ub'+str(i)+'.top -o '+args.wd+'/ub'+str(i)+'_em.tpr')
    if args.em:
        # for i in 1, 2, 4, 8:
        i=1
        os.system(args.gromacs_command+
            ' mdrun -deffnm '+args.wd+'/ub'+str(i)+'_em')
        for rep in 1, 2, 3:
            os.system(args.gromacs_command+
                ' grompp -f MDP/nvt'+mdp_index+'.mdp -c '+args.wd+'/ub'+str(i)+'_em.gro -p '+args.wd+'/ub'+str(i)+'.top -o '+args.wd+'/ub'+str(i)+'_nvt.'+str(rep)+'.tpr')
    if args.nvt:
        for i in 1, 2, 4, 8:        
            for rep in 1, 2, 3:
                # os.system(args.gromacs_command+
                #     ' mdrun -deffnm '+args.wd+'/ub'+str(i)+'_nvt.'+str(rep))
                os.system(args.gromacs_command+
                    ' grompp -f MDP/npt'+mdp_index+'.mdp -c '+args.wd+'/ub'+str(i)+'_nvt.'+str(rep)+'.gro -p '+args.wd+'/ub'+str(i)+'.top -o '+args.wd+'/ub'+str(i)+'_npt.'+str(rep)+'.tpr')
    if args.npt:
        for i in 1, 2, 4, 8:
            for rep in 1, 2, 3:
                # os.system(args.gromacs_command+
                #     ' mdrun -deffnm '+args.wd+'/ub'+str(i)+'_npt.'+str(rep))
                os.system(args.gromacs_command+
                    ' grompp -f MDP/prod'+mdp_index+'.mdp -c '+args.wd+'/ub'+str(i)+'_npt.'+str(rep)+'.gro -p '+args.wd+'/ub'+str(i)+'.top -o '+args.wd+'/ub'+str(i)+'_prod.'+str(rep)+'.tpr')    

    if args.analyze:
        for i in 1, 2, 4, 8:
            os.system('cp '+args.wd+'/ub'+str(i)+'_npt.gro '+args.wd+'/ub'+str(i)+'_npt.3.gro')
            for rep in 1, 2, 3:            
                os.system(args.gromacs_command+
                    ' grompp -f MDP/analyze'+mdp_index+'.mdp -c '+args.wd+'/ub'+str(i)+'_npt.'+str(rep)+'.gro -p '+args.wd+'/ub'+str(i)+'.top -o '+args.wd+'/ub'+str(i)+'_analyze.'+str(rep)+'.tpr')    
            

 #    if args.disre:
 #        # for i in 1, 2, 4, 8:
 #        i=1
 #        os.system(args.gromacs_command+
 #            ' grompp -f MDP/disre'+mdp_index+'.mdp -c '+args.wd+'/ub'+str(i)+'_prod_protein.gro -p '+args.wd+'/ub'+str(i)+'_disre.top -o '+args.wd+'/ub'+str(i)+'_disre.tpr')

 # 