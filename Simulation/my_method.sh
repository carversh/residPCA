module load Anaconda2/5.2.0
source activate R_4.1.0_python_3.8.2
# module load Anaconda3/5.2.0
# source activate R_4.1.0_python_3.8.2
#LSB_JOBINDEX=$1
# load in the parameters
param_line=`cat my_method.par | head -n $LSB_JOBINDEX | tail -n1`
param_line=`cat across_cts_40000.par | head -n $LSB_JOBINDEX | tail -n1`

PARAM1=`echo $param_line | awk '{ print $1 }'`
PARAM2=`echo $param_line | awk '{ print $2 }'`
PARAM3=`echo $param_line | awk '{ print $3 }'`
PARAM4=`echo $param_line | awk '{ print $4 }'`
PARAM5=`echo $param_line | awk '{ print $5 }'`
PARAM6=`echo $param_line | awk '{ print $6 }'`
PARAM7=`echo $param_line | awk '{ print $7 }'`
PARAM8=`echo $param_line | awk '{ print $8 }'`

echo $PARAM1 $PARAM2 $PARAM3 $PARAM4 $PARAM5 $PARAM6 $PARAM7 $PARAM8

Rscript Simulation.r $PARAM1 $PARAM2 $PARAM3 $PARAM4 $PARAM5 $PARAM6 $PARAM7 $PARAM8

conda deactivate
source activate scanpy_env_eris2

python my_method.py $PARAM1 $PARAM2 $PARAM3 $PARAM4 $PARAM5 $PARAM6 $PARAM7 --flag_sim_100_cts $PARAM8

















#python my_method.py 1 my_method within_one_ct 105 40000 0.15 0.15 --flag_sim_100_cts cts_100
#echo $PARAM1 $PARAM2 $PARAM3 $PARAM4 $PARAM5 $PARAM6 $PARAM7 $PARAM8 >> finished_jobs_20k_.txt

# bsub -o OUT.%J.%I.out -q bigmem -R 'rusage[mem=33000]' -W 2:00 -J ct100[962-1281] "bash simulation_erisTWO.sh"
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=35000]' -W 6:00 -J ct100[321-640] "bash my_method.sh"

# bsub -o OUT.%J.%I.out -q rerunnable -R 'rusage[mem=12000]' -W 2:00 -J simulation[1-320] "bash my_method.sh"
# bsub -o OUT.%J.%I.out -q rerunnable -R 'rusage[mem=12000]' -W 5:00 -J simulation[321-640] "bash simulation_erisTWO.sh"

# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=50000]' -W 60:00 -J nmfsim[1921] "bash simulation_erisTWO.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=50000]' -W 60:00 -J nmfsim[2241] "bash simulation_erisTWO.sh"

# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=30000]' -W 85:00 -J nmfsim[224-320] "bash simulation_erisTWO.sh"



