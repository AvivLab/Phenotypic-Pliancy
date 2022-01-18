#!/bin/bash
### shell type
#$ -S /bin/bash
### current working directory as the starting place
#$ -cwd
### name of standard output  file
#$ -o dataCollection_out
### name of standar error file
#$ -e dataCollection_err
### whether to join error and output files (y,n)
#$ -j n
###
#$ -r y
### name of project
#$ -N collection
#$ -q highmem.q
#$ -l h_vmem=20G
### number of tasks
#
############# $ -t 1-310360


# using local julia installation in home folder julia 1.6.4
cd /gs/gsfs0/users/avibergm/Polycomb/PolycombCode/PolycombPaper/src.1
CMDFILE=/gs/gsfs0/users/avibergm/Polycomb/PolycombCode/PolycombPaper/src.1/dataCollectionInput.txt
CMD=$(awk "NR==$SGE_TASK_ID" $CMDFILE)
$CMD
