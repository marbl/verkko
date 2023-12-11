import os
import verkkoHelper as vH
import snakemake.io as sIO


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


def run(inputs, outputs, log, params, threads, resources, wildcards):
  os.makedirs('1-prepare', exist_ok=True)

  with open(outputs.out1, "w") as f1:
    for i1 in range(10,20):
      f1.write(f'{i1}\n')

  with open(outputs.out2, "w") as f2:
    for i2 in range(20,30):
      print(f'{i2}', file=f2)

