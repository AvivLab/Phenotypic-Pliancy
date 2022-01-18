#!/bin/bash
### shell type
#$ -S /bin/bash
### current working directory as the starting place
#$ -cwd
### name of standard output  file
#$ -o polycomb_out
### name of standar error file
#$ -e polycomb_err
### whether to join error and output files (y,n)
#$ -j n
###
#$ -r y
### name of project
#$ -N polycomb
### number of tasks
#
############# $ -t 1-310360

# load julia module
#module load julia/1.0.4
#module load julia/1.2.0 # lets use local julia installation (home folder)
#cd /home/aviv/Polycomb/julia1.1.1_MarylNewVersionWithTwoEnvs_Aug2019/src/

module load julia/1.1.1

# using local julia installation in home folder julia 1.1.1 
cd /gs/gsfs0/users/avibergm/Polycomb/PolycombCode/PolycombPaper/src
CMDFILE=/gs/gsfs0/users/avibergm/Polycomb/PolycombCode/PolycombPaper/src/parameters.txt
CMD=$(awk "NR==$SGE_TASK_ID" $CMDFILE)
$CMD





