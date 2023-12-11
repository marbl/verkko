import os
import verkkoHelper as vH
import snakemake.io as sIO

import verkkoA


def inputs(rules, wildcards, checkpoints):
  outA = vH.ruleOutputs(verkkoA)
  return { 'in1': outA.out1,
           'in2': outA.out2 }

def outputs(rules):
  return { 'finished': '2-partition/finished' }

def logs(rules):
  return {}

def params(rules):
  return {}

def threads(rules):
  return 1

def resources(rules):
  return {}


def run(inputs, outputs, log, params, threads, resources, wildcards):
  os.makedirs('2-partition', exist_ok=True)
  os.makedirs('3-compute', exist_ok=True)

  f1  = open(inputs.in1)
  f2  = open(inputs.in2)
  idx = 0

  l1 = f1.readline()
  l2 = f2.readline()
  while l1 and l2:
    idx = idx+1

    ot = open(f'2-partition/part{format(idx, "03d")}.txt', 'w')
    print(f'{l1}{l2}', file=ot)
    ot.close()

    l1 = f1.readline()
    l2 = f2.readline()
 
  f1.close()
  f2.close()

  fin = open(outputs.finished, "w")
  fin.close()
