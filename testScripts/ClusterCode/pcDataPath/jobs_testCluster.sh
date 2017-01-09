mcc -m /home/mpib/perdikis/testClusterFolder/testCluster.m -a /home/mpib/perdikis/testClusterFolder/Toolboxes
for iJ in {1..3}; do
echo "#PBS -N test$iJ" >> jobfile
echo "#PBS -l walltime=24:00:00" >> jobfile
echo "#PBS -l mem=1gb" >> jobfile
echo "#PBS -j oe" >> jobfile
echo "/home/mpib/perdikis/testClusterFolder/run_testCluster.sh /opt/matlab/interactive $iJ" >> jobfile
qsub jobfile
rm -f jobfile
done