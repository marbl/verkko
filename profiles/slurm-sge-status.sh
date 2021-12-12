#!/bin/sh

#  The status.sh script must return the status of the supplied grid job id.
#    running, success, or failed.

#  JobID and return status.
jobid=$1
jobstatus="failed"

#  Test for Slurm: if sinfo is present, assume slurm works.
slurm=`which sinfo 2> /dev/null`

#  Test for SGE: if SGE_CELL exists in the environment, assume SGE works.
sge=$SGE_CELL


#  Check Slurm status.
#
#  Derived from https://github.com/jdblischak/smk-simple-slurm
#
#  Note that if sacct is queried immediately after job submission, it
#  is possible that sacct returns no output, and so our default state
#  is "running".
#
if [ "x$slurm" != "x" ] ; then
  jobstatus=`sacct -j "$jobid" --format State --noheader | head -n 1 \
  | \
  awk \
  'BEGIN { stat="running" } \
   { \
    if      ($1 == "COMPLETED")   { stat="success" } \
    else if ($1 == "RUNNING")     { stat="running" } \
    else if ($1 == "PENDING")     { stat="running" } \
    else if ($1 == "COMPLETING")  { stat="running" } \
    else if ($1 == "CONFIGURING") { stat="running" } \
    else if ($1 == "SUSPENDED")   { stat="running" } \
    else                          { stat="failed"  }
   }
   END { print stat }'`


#  Check SGE status.
#
#  Query qstat for the job.  If the job is in an error state, it will report
#  'error reason', otherwise, if there is output, the job is pending or
#  running.
#
#  If no output, then we need to query qacct for the job status.
#
elif [ "x$sge" != "x" ] ; then
  if   [ `qstat -j $jobid 2> /dev/null | grep error\ reason | wc -c` -gt 0 ] ; then
    jobstatus="failed"

  elif [ `qstat -j $jobid 2> /dev/null | wc -c` -gt 0 ] ; then
    jobstatus="running"

  else
    jobstatus=`qacct -j $jobid 2> /dev/null \
    | \
    awk \
    '{ \
      if (($1 == "failed")      && ($2 > 0)) { fail=$2 } \
      if (($1 == "exit_status") && ($2 > 0)) { fail=$2 } \
     } END { \
      if (fail > 0) { print "failed" }  \
      else          { print "success" } \
     }'`
  fi

#  Otherwise, do what?  Fail!
else
    jobstatus="failed"
fi

echo $jobstatus
exit 0
