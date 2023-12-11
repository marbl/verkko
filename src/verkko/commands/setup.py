import os                    #
import sys                   #
import inspect               #  For 'cleandoc'
#mport configparser as cP    #  For reading verkko.ini
import vconfig
import verkkoHelper as vH    #  For 'ruleOutput' and 'ruleInputs'
import snakemake.io as sIO   #  For 'glob_wildcards' and 'expand'

#
#  See comments in verkko.py, please.
def synopsis():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''validate command line options and write a config file for the assembly'''
  return inspect.cleandoc(s)

def details():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''The 'setup' command creates a parameter file ('verkko.ini') in the run-directory
        which contains all parameters and input file names for the assembly.

        It will check that:
         - parameters are valid
         - parameters are consistent with any existing parameter file
         - input files exist and are readable
         - external programs exist and can be executed

        Any failing check will result in an immediate failure, except that
        when parameters are not consistent with the existing parameter file,
        it will:
         - silently update the parameters if no intermediate results depend
           on the changed parameters
         - otherwise, depending on option '--recompute', emit an error indicating
           which parameter changed and which intermediate result(s) are affected,
           or remove all such intermediate results and recomptue.

        All other Verkko sub-commands read their parameters from 'verkko.ini'.
        Suuplying algorithmic parameters directly to the sub-command will override
        the global parameter but NOT update the parameter file.'''
  return inspect.cleandoc(s)



#def inputs(rules, wildcards, checkpoints):
#  return {}
#
#def outputs(rules):
#  return { 'ini': 'verkko.current.ini' }
#
#def logs(rules):
#  return {}
#
#def params(rules):
#  return {}
#
#def threads(rules):
#  return 1
#
#def resources(rules):
#  return {}
#


def run():
  cfg = vconfig.verkkoConfig()
  cfg.load('verkko.current.ini')

  checkValid(cfg)
  checkInputs(cfg)
  checkExternal(cfg)
  checkConfig(cfg)
  writeConfig(cfg)


class cCint(object):
  def __init__(self, val):              self.val = int(val)
  def __enter__(self):                  return self.val
  def __exit__(self, type, value, tb):  pass

class cCfloat(object):
  def __init__(self, val):              self.val = float(val)
  def __enter__(self):                  return self.val
  def __exit__(self, type, value, tb):  pass

def fileExistsAndIsReadable(src, fn):
  if not os.path.exists(fn):
    return(f"input {src} '{fn}' does not exist\n")
  return(f'')




def checkValid(cfg, cErrors):

  with cCint(cfg['CANU_READ_CORRECTION']['kmer-size']) as p:
    if p < 14 or p > 32:
      cErrors.append(f'Invalid --correct-kmer-size {p}: must be between 14 and 32 inclusive')

  #with cCint(cfg['CANU_READ_CORRECTION']['kmer-threshold']) as p:
  #  if p < 14 or p > 32:
  #    cErrors.append(f'Invalid --correct-kmer-size {p}: must be between 14 and 32 inclusive')

  #with cCfloat(cfg['CANU_READ_CORRECTION']['kmer-filter']) as p:
  #  if p < 14 or p > 32:
  #    cErrors.append(f'Invalid --correct-kmer-size {p}: must be between 14 and 32 inclusive')

  #with cCint(cfg['CANU_READ_CORRECTION']['min-read-length']) as p:
  #  if p < 14 or p > 32:
  #    cErrors.append(f'Invalid --correct-kmer-size {p}: must be between 14 and 32 inclusive')

  #with cCint(cfg['CANU_READ_CORRECTION']['min-overlap-length']) as p:
  #  if p < 14 or p > 32:
  #    cErrors.append(f'Invalid --correct-kmer-size {p}: must be between 14 and 32 inclusive')

  #with cCint(cfg['CANU_READ_CORRECTION']['min-kmers']) as p:
  #  if p < 14 or p > 32:
  #    cErrors.append(f'Invalid --correct-kmer-size {p}: must be between 14 and 32 inclusive')

  #with cCint(cfg['CANU_READ_CORRECTION']['hash-bits']) as p:
  #  if p < 14 or p > 32:
  #    cErrors.append(f'Invalid --correct-kmer-size {p}: must be between 14 and 32 inclusive')


def checkInputs(cfg, cErrors):
  errors = ''

  for fn in cfg['HIFI'].values:
    print(f'DumpChecking HiFi:  \'{fn}\'')
  for fn in cfg['ONT'].values:
    print(f'DumpChecking ONT:   \'{fn}\'')

  for fn in cfg['HIFI'].values:
    print(f'Checking HiFi:  \'{fn}\'')
    errors += fileExistsAndIsReadable('HiFi', fn)

  for fn in cfg['ONT'].values:
    print(f'Checking ONT:   \'{fn}\'')
    errors += fileExistsAndIsReadable('ONT', fn)

  for fn in cfg['HIC1'].values:
    print(f'Checking HiC-2: \'{fn}\'')
    errors += fileExistsAndIsReadable('HiC-1', fn)

  for fn in cfg['HIC2'].values:
    print(f'Checking HiC-1: \'{fn}\'')
    errors += fileExistsAndIsReadable('HiC-2', fn)

  #for fn in cfg['HAP1']:
  #  errors += fileExistsAndIsReadable(fn)

  #for fn in cfg['HAP2']:
  #  errors += fileExistsAndIsReadable(fn)

  print(errors)

def checkExternal(cfg, cErrors):
  a=3

def checkConfig(cfg, cErrors):
  a=3

def writeConfig(cfg, cErrors):
  a=3


#needs a config step that
# - compares existing options vs config files
# - writes the list of input reads to config files
# - writes all parameters to config files
# - later steps read ALL parameters/inputs from these files, NOT COMMAND LINE
