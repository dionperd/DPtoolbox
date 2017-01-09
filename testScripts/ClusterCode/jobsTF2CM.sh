#!/bin/bash

mcc -m KHVATF2CM_grid.m -a Toolboxes/

for ct in {1..3}; do

#set the names of the jobs
		var="CT$ct"
		echo "#PBS -N TF2CM-CT$var" >> jobfile
#set the resources of the jobs
		echo "#PBS -l walltime=30:00:00" >> jobfile
		echo "#PBS -l mem=50gb" >> jobfile
#set the names of the jobs
		echo "#PBS -j oe" >> jobfile
		echo "./run_KHVATF2CM_grid.sh /opt/matlab/interactive $ct" >> jobfile
		qsub jobfile
		rm -f jobfile
done