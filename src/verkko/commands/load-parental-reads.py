import sys
import inspect

#
#  Sub-command description for user-interface.
#    (see comments in verkko.py, please.)
#

def synopsis():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''load parental reads into the Verkko run directory'''
  return inspect.cleandoc(s)

def details():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''Sub-command load-parental-reads copies parental Illumina reads
        from their source location to the Verkko assembly directory.'''
  return inspect.cleandoc(s)

#
#  Workflow description and parameters.
#

def inputs(rules, wildcards, checkpoints):
  return {}

def outputs(rules):
  return { 'out1': '1-prepare/Aseq1',
           'out2': '1-prepare/Aseq2' }

def logs(rules):
  return {}

def params(rules):
  return {}

def threads(rules):
  return 1

def resources(rules):
  return {}

#
#  Sub-command execution.
#

def run(args):
  print(args)

