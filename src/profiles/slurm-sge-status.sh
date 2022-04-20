#!/bin/sh

#  The status.sh script must return the status of the supplied grid job id.
#    running, success, or failed.

#  JobID and return status.
jobid=$1
jobstatus="running"

#  Test for Slurm: if sinfo is present, assume slurm works.
slurm=$(which sinfo 2> /dev/null)

#  Test for SGE: if SGE_ROOT exists in the environment, assume SGE works.
sge=$SGE_ROOT

#  Test for LSF: if LSF_ENVDIR exists in the environment, assume LSF works.
lsf=$LSF_ENVDIR


#  Check Slurm status.
#
#  Derived from https://github.com/jdblischak/smk-simple-slurm
#
#  Note that if sacct is queried immediately after job submission, it
#  is possible that sacct returns no output, and so our default state
#  is "running".
#
#  dashboard_cli is an NIH biowulf extension that caches sacct results.
#
if [ "x$slurm" != "x" ] ; then
  if [ -x /usr/local/bin/dashboard_cli ] ; then
    sleep 2
    sacct="/usr/local/bin/dashboard_cli jobs --fields state"
  else
    sacct="sacct --format State"
  fi

  jobstatus=$($sacct -j "$jobid" --noheader)
  if [ $? != 0 ] ; then   #  If sacct fails, report
    echo running          #  that the job is running.
    echo 1>&2 "Job $jobid is $jobstatus (by default)."
    exit 0
  fi

  jobstatus=$(echo $jobstatus | head -n 1 \
  | \
  awk \
   'BEGIN { stat="running" } \
   { \
    if      ($1 == "COMPLETED")     { stat="success" } \
    else if ($1 == "PENDING")       { stat="running" } \
    else if ($1 == "CONFIGURING")   { stat="running" } \
    else if ($1 == "RUNNING")       { stat="running" } \
    else if ($1 == "SUSPENDED")     { stat="running" } \
    else if ($1 == "PREEMPTED")     { stat="running" } \
    else if ($1 == "COMPLETING")    { stat="running" } \
    else if ($1 == "BOOT_FAIL")     { stat="failed"  } \
    else if ($1 == "CANCELLED")     { stat="failed"  } \
    else if ($1 == "DEADLINE")      { stat="failed"  } \
    else if ($1 == "FAILED")        { stat="failed"  } \
    else if ($1 == "NODE_FAIL")     { stat="failed"  } \
    else if ($1 == "OUT_OF_MEMORY") { stat="failed"  } \
    else if ($1 == "PREEMPTED")     { stat="failed"  } \
    else if ($1 == "TIMEOUT")       { stat="failed"  } \
    else                            { stat="failed"  }
   }
   END { print stat }')


#  Check SGE status.
#
#  Query qstat for the job.  If the job is in an error state, it will report
#  'error reason', otherwise, if there is output, the job is pending or
#  running.
#
#  If no output, then we need to query qacct for the job status.
#
elif [ "x$sge" != "x" ] ; then
  if   [ $(qstat -j $jobid 2> /dev/null | grep error\ reason | wc -c) -gt 0 ] ; then
    jobstatus="failed"

  elif [ $(qstat -j $jobid 2> /dev/null | wc -c) -gt 0 ] ; then
    jobstatus="running"

  else
    jobstatus=$(qacct -j $jobid 2> /dev/null \
    | \
    awk \
    '{ \
      if (($1 == "failed")      && ($2 > 0)) { fail=$2 } \
      if (($1 == "exit_status") && ($2 > 0)) { fail=$2 } \
     } END { \
      if (fail > 0) { print "failed" }  \
      else          { print "success" } \
     }')
  fi


#  Check LSF status.
#
#  https://github.com/Snakemake-Profiles/lsf/blob/master/%7B%7Bcookiecutter.profile_name%7D%7D/lsf_status.py
#    running -- PEND RUN PSUSP USUSP SSUSP WAIT
#    success -- DONE POST_DONE
#    failed  -- UNKWN ZOMBI EXIT POST_ERR
#
elif [ "x$lsf" != "x" ] ; then
  jobstatus=$(bjobs -o stat -noheader "$jobid" | head -n 1 \
  | \
  awk \
  'BEGIN { stat="running" } \
   { \
    if      ($1 == "DONE")        { stat="success" } \
    else if ($1 == "POST_DONE")   { stat="success" } \
    else if ($1 == "PEND")        { stat="running" } \
    else if ($1 == "RUN")         { stat="running" } \
    else if ($1 == "PSUSP")       { stat="running" } \
    else if ($1 == "USUSP")       { stat="running" } \
    else if ($1 == "SSUSP")       { stat="running" } \
    else if ($1 == "WAIT")        { stat="running" } \
    else if ($1 == "UNKWN")       { stat="failed"  } \
    else if ($1 == "ZOMBI")       { stat="failed"  } \
    else if ($1 == "EXIT")        { stat="failed"  } \
    else if ($1 == "POST_ERR")    { stat="failed"  } \
    else                          { stat="failed"  } \
   } \
   END { print stat }')


#  Otherwise, do what?  Fail!
else
    jobstatus="failed"
fi

#echo 1>&2 "Job $jobid is $jobstatus."

echo $jobstatus
exit 0
