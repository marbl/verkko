import os
import verkkoHelper as vH
import snakemake.io as sIO

import verkkoB


def inputs(rules, wildcards, checkpoints):
  outB = vH.ruleOutputs(verkkoB)
  return { 'par': outB.finished,
           'inp': '2-partition/part{xxxx}.txt' }

def outputs(rules):
  return { 'out': '3-compute/part{xxxx}.xtx' }

def logs(rules):
  return {}

def params(rules):
  return {}

def threads(rules):
  return 2

def resources(rules):
  return {}


def run(inputs, outputs, log, params, threads, resources, wildcards):

  inp = open(inputs.inp)
  out = open(outputs.out, 'w')

  for f in inp:
    print(f'{f} {f}', file=out)
