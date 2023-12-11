import os
#rom   snakemake.io import glob_wildcards, expand
from   collections  import namedtuple

def ruleOutputs(mod):
  return namedtuple('X', mod.outputs({}).keys())(*mod.outputs({}).values())

def ruleInputs(mod):
  return namedtuple('X', mod.inputs({}, {}, {}).keys())(*mod.inputs({}, {}, {}).values())
