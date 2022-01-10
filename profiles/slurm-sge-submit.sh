#!/bin/sh

#  This script is responsible for submitting jobs to the grid.  The job
#  will be the last argument; other args will be as specified in
#  the --cluster option to snakemake.
#
#  The sript must print the job ID to stdout and only that.
#
#  No rewriting of this script is performed by snakemake (unlike the
#  --jobscript template), but the --cluster option can pass in resources to
#  use:
#    -cluster "./cluster.sh {threads} {resources.mem_gb} {resources.anything}"
#  will pass in the 'threads' value from the rule, and 'mem_gb' and 'anything'
#  entries in resources.  ALL rules must have ALL resources specified.
#
#  Default resources can be specified using '--default-resources', or
#  as an array in the profile config file:
#    default-resources:
#     - blah=200
#     - mem_gb=100
#

#  We expect (as defined in the --cluster option) the first three args to be
#  as below, then any number of user-supplied args, then the name of the
#  script to run.
#
#  Since we want the user-supplied args to be at the end of the sbatch/qsub
#  arg list, and since the script to run is expected to be the last arg to
#  sbatch/qsub, we can extract the three values below and then simply pass
#  the rest of our args to sbatch/qsub.

n_cpus=$1; shift   #  Number of CPUs we can use
mem_gb=$1; shift   #  Exepcted memory needed, in gigabytes
time_h=$1; shift   #  Expected time needed, in hours
rule_n=$1; shift   #  Name of the rule
jobidx=$1; shift   #  Index of the job, "1" if a non-parallel job

eval script=\${$#} #  Path to the script we want to submit

#  Test for Slurm: if sinfo is present, assume slurm works.
slurm=`which sinfo 2> /dev/null`

#  Test for SGE: if SGE_CELL exists in the environment, assume SGE works.
sge=$SGE_CELL

#  Test for LSF: if LSF_ENVDIR exists in the environment, assume LSF works.
lsf=$LSF_ENVDIR


#  Submit to Slurm.
#
#  Note that --time expects format hh:mm:ss.
#
if [ "x$slurm" != "x" ] ; then
  jobid=$(sbatch --cpus-per-task ${n_cpus} --mem ${mem_gb}g --time ${time_h}:00:00 --output batch-scripts/%A.${rule_n}.${jobidx}.out "$@")


#  Submit to SGE.  '-terse' reports just the job ID instead of the usual
#  human-parsable text.
#
#  Note that SGE expects memory as "memory per thread" instead of
#  "memory per task" that we're supplied.
#
elif [ "x$sge" != "x" ] ; then
  mem_per_thread=$(dc -e "3 k ${mem_gb} ${n_cpus} / p")
  jobid=$(qsub -terse -cwd -V -pe thread ${n_cpus} -l memory=${mem_per_thread}g -j y -o batch-scripts/\$JOB_ID.${rule_n}.${jobidx}.out "$@")


#  Submit to LSF.
#  Other options:
#    -A account
#    -q queue
#    -p partition
#    -n ?
#    -J job name
#    -W time limit
#    -R resources needed -- "\"select[mem>20000] rusage[mem=20000] span[hosts=1]\"",
#    -e err-out
#    -o out-out
#
elif [ "x$lsf" != "x" ] ; then
  mem_units=$(grep ^LSF_UNIT_FOR_LIMITS $LSF_ENVDIR/lsf.conf | cut -d= -f 2)

  if   [ $mem_units = "GB" ] ; then
    mem=${mem_gb}
  elif [ $mem_units = "MB" ] ; then
    mem=$(dc -e "3 k ${mem_gb} 1024 / p")
  else
    mem=$(dc -e "3 k ${mem_gb} 1048576 / p")
  fi

  jobid=$(bsub -eo batch-scripts/\$JOB_ID.${rule_n}.${jobidx}.err \
               -oo batch-scripts/\$JOB_ID.${rule_n}.${jobidx}.out \
               -R span[hosts=1] -n ${n_cpus} -M ${mem} "$@")


#  Otherwise, fail.
else
  exit 1
fi

#  Save a copy of the script we run, for debugging.
mkdir -p        batch-scripts/
cp -p ${script} batch-scripts/${jobid}.${rule_n}.${jobidx}.sh

#  Report the job ID back to snakemake.
echo $jobid
exit 0
