#!/bin/sh

#  This script forms a template for the job that is eventually run.
#  {{exec_job}} at the end will become a call to snakemake to run the desired
#  job.
#
#  Some things you can get here:
#    {{properties}}  - a text string of all the properties for this job
#    {{job}}         - the name of the rule running
#    {{jobfinished}} - name of the sentinel file announcing a successful job
#    {{jobfailed}}   - name of the sentinel file announcing a failure
#    {{target}}      - name of the output file(s)
#    {{cores}}       - 'all'
#
#    {{workflow.workdir_init}}   - os.path.abspath(os.curdir)
#    {{workflow.attempt}}        - number of attempts to run this job
#    {{workflow.jobscript}}
#
#    {{workflow.cores}}          - number of cores to use
#    {{workflow.nodes}}          - number of nodes to use
#
#    {{workflow.basedir}}        - os.path.dirname(snakefile)
#    {{workflow.main_snakefile}} - os.path.abspath(snakefile)
#

#  Save a copy of the script we run, for debugging.
#mkdir -p {workflow.workdir_init}/scripts/
#cp -p $0 {workflow.workdir_init}/scripts/

#  Run the job.
{exec_job}
