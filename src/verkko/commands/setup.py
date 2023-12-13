import os                    #
import sys                   #
#mport inspect               #  For 'cleandoc'
#mport configparser as cP    #  For reading verkko.ini
import verkkoConfig as vC    #
#mport verkkoHelper as vH    #  For 'ruleOutput' and 'ruleInputs'
#mport snakemake.io as sIO   #  For 'glob_wildcards' and 'expand'

#
#  Validate the parameter file ('verkko.current.ini') in the run-directory
#  which contains all parameters and input file names for the assembly.
#
#  It will check that:
#   - parameters are valid
#   - parameters are consistent with any existing parameter file
#   - input files exist and are readable
#   - external programs exist and can be executed
#
#  Any failing check will result in an immediate failure, except that
#  when parameters are not consistent with the existing parameter file,
#  it will:
#   - silently update the parameters if no intermediate results depend
#     on the changed parameters
#   - otherwise, depending on option '--recompute', emit an error indicating
#     which parameter changed and which intermediate result(s) are affected,
#     or remove all such intermediate results and recomptue.
#
#  All other Verkko sub-commands read their parameters from 'verkko.ini'.
#  Suuplying algorithmic parameters directly to the sub-command will override
#  the global parameter but NOT update the parameter file.'''
#


#  Helpers.  Convert a string to an integer/float and
#  allow use in a 'with cCint(string) as p' block.
#
class cCint(object):
  def __init__(self, val):              self.val = int(val)
  def __enter__(self):                  return self.val
  def __exit__(self, type, value, tb):  pass

class cCfloat(object):
  def __init__(self, val):              self.val = float(val)
  def __enter__(self):                  return self.val
  def __exit__(self, type, value, tb):  pass


#  Test if a file exists, and if not, return an error string.
def fileExistsAndIsReadable(src, fn, err):
  if not os.path.exists(fn):
    err.append(f'{src} input \'{fn}\' does not exist')



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
  for fn in cfg['HIFI'].values:   fileExistsAndIsReadable('HiFi',  fn, cErrors)
  for fn in cfg['ONT' ].values:   fileExistsAndIsReadable('ONT',   fn, cErrors)
  for fn in cfg['HIC1'].values:   fileExistsAndIsReadable('HiC-1', fn, cErrors)
  for fn in cfg['HIC2'].values:   fileExistsAndIsReadable('HiC-2', fn, cErrors)
  for fn in cfg['HAP1'].values:   fileExistsAndIsReadable('Hap-1', fn, cErrors)
  for fn in cfg['HAP2'].values:   fileExistsAndIsReadable('Hap-2', fn, cErrors)

def checkExternal(cfg, cErrors):
  a=3

def checkConfig(cfg, cErrors):
  a=3

def failIfErrors(cfg, cErrors):
  if len(cErrors) == 0:
    return

  print(f'')
  print(f'Errors detected:')
  for e in cErrors:
    print(f' - {e}')
  print(f'')
  exit(1)

#needs a config step that
# - compares existing options vs config files
# - writes the list of input reads to config files
# - writes all parameters to config files
# - later steps read ALL parameters/inputs from these files, NOT COMMAND LINE
