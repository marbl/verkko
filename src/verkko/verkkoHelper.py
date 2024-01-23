import os
#rom   snakemake.io import glob_wildcards, expand
from   collections  import namedtuple
import subprocess

#
#  SnakeMake stuff.
#
#  ruleOutputs() returns a dict of the output for the specified rule:
#       x = ruleOutput('importReads')
#       doSomething(x.reads)
#
#  ruleInputs() is similar, but for the inputs to a rule.
#

def ruleOutputs(mod):
  return namedtuple('X', mod.outputs({}).keys())(*mod.outputs({}).values())

def ruleInputs(mod):
  return namedtuple('X', mod.inputs({}, {}, {}).keys())(*mod.inputs({}, {}, {}).values())

#
#  A nicer exception.  It gets caught at the end of verkko.py, catching
#  (nearly) all exceptions thrown.  If you raise a verkkoException and pass
#  in a list of strings (or a list of lists, etc), they're printed
#  one-per-line.
#
class verkkoException(Exception):
  def __init__(self, messages):
    self.message = '\n  '.join(flatten(messages))

  def __str__(self):
    return self.message


#
#  Directory management.
#    enterdir(dirname)
#     - save the current directory in 'dirstack'.
#     - create, if needed, 'dirname' and chdir() into it.
#
#    leavedir()
#     - return to whatever directory was last pushed onto 'dirstack'.
#
dirstack = []

def createdir(dir):
  try:
    os.makedirs(dir, exist_ok=True)
  except OSError as e:
    od = os.getcwd()
    raise verkkoException([f'In directory \'{od}\':',
                           f'Cannot create directory \'{e.filename}\': {e.strerror}'])

def enterdir(dir):
  od = os.getcwd()
  nd = od + '/' + dir
  createdir(nd)
  try:
    dirstack.append(od)
    createdir(dir)
    os.chdir(dir)
    nd = os.getcwd()
    print(f'Entering \'{nd}\' from \'{od}\'')
  except OSError as e:
    raise verkkoException([f'In directory \'{od}\':',
                           f'Cannot use directory \'{e.filename}\': {e.strerror}'])

def leavedir():
  od = os.getcwd()
  nd = 'nowhere'
  if len(dirstack) == 0:
    raise verkkoException(f'Attempt to leavedir() from \'{od}\' to <nowhere>')
  try:
    nd = dirstack.pop()
    os.chdir(nd)
    print(f'Leavng \'{od}\' returning to \'{nd}\'')
  except OSError as e:
    raise verkkoException([f'ERROR: Cannot return to  directory -d \'{e.filename}\': {e.strerror}'])

#
#  Recursively convert a thing (a scalar, list or tuple) to a list.
#

def flattenlist(items):
  for item in items:
    if isinstance(item, list) or isinstance(item, tuple):
      yield from flattenlist(item)
    else:
      yield item

def flatten(thing):
  output = []
  if isinstance(thing, list) or isinstance(thing, tuple):
    output += flattenlist(thing)
  else:
    output.append(thing)
  return output

