import os
import verkkoHelper as vH
import snakemake.io as sIO

import verkkoA
import verkkoB
import verkkoC


def inputs(rules, wildcards, checkpoints):
  outA = vH.ruleOutputs(verkkoA)
  outB = vH.ruleOutputs(verkkoB)
  inpC = vH.ruleInputs (verkkoC)
  outC = vH.ruleOutputs(verkkoC)

  print(f'inpC {inpC}')
  print(f'outC {outC}')

  w = sIO.glob_wildcards(inpC.inp).xxxx     #  Find wildcard matches to input files.
  p = sIO.expand        (outC.out, xxxx=w)  #  Replace wildcards in outputs.

  return { 'in1':         outA.out1,
           'in2':         outA.out2,
           'partitioned': outB.finished,
           'pieces':      p }

def outputs(rules):
  return { 'merged': '4-output/merged' }

def logs(rules):
  return {}

def params(rules):
  return {}

def threads(rules):
  return 1

def resources(rules):
  return {}


def run(inputs, outputs, log, params, threads, resources, wildcards):
  os.makedirs('4-output', exist_ok=True)

  out = open(outputs.merged, 'w')
  for f in inputs.pieces:
    print(f'piece {f}', file=out)
  out.close()
