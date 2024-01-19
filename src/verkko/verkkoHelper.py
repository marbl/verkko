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
#  (nearly) all exceptions thrown.  If you raise a verkkoException and
#  pass in a list of strings, they're printed one-per-line.
#
#
class verkkoException(Exception):
  def __init__(self, lines):
    self.message = ''
    if isinstance(lines, str):
      self.message = lines
    else:
      self.message = '\n  '.join(lines)

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

def enterdir(dir):
  od = os.getcwd()
  nd = od + '/' + dir
  x = 0/0
  try:
    dirstack.append(od)
    os.makedirs(dir, exist_ok=True)
    os.chdir(dir)
    nd = os.getcwd()
    print(f'Entering \'{nd}\' from \'{od}\'')
  except FileExistsError:
    raise verkkoException([f'In directory \'{od}\':',
                           f'Cannot create directory -d \'{e.filename}\': {e.strerror}'])
  except OSError as e:
    raise verkkoException([f'In directory \'{od}\':',
                           f'Cannot create or use directory -d \'{e.filename}\': {e.strerror}'])

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

#
#  Execute an external command.
#   - command is an absolute path to the binary to run
#   - args is a list of parameters to pass
#   - stdin  must be from an actual file (or /dev/null)
#   - stdout must be  to  an actual file (or /dev/null)
#   - stderr must be  to  an actual file (or /dev/null)
#   - no shell interpretation is performed
#      - pipelines are not supported
#      - file name globbing not supported
#
#  Input/output from/to a PIPE is NOT ENCOURAGED.  It is exceptionally easy
#  to deadlock.  Deadlock is all but guaranteed if more than one PIPE is used.
#

def runExternal(command, *args, stdin=None, stdout=None, stderr=None):
  if stdin  == None or stdin  == '/dev/null':   stdin  = subprocess.DEVNULL
  if stdout == None or stdout == '/dev/null':   stdout = subprocess.DEVNULL
  if stderr == None or stderr == '/dev/null':   stderr = subprocess.DEVNULL

  if isinstance(stdin,  str):    stdin = open(stdin,  'r')
  if isinstance(stdout, str):   stdout = open(stdout, 'w')
  if isinstance(stderr, str):   stderr = open(stderr, 'w')

  cmdarg = [command] + flatten(args)
  cmddis = cmdarg.copy()
  cmddis.reverse()

  print(f'----------------------------------------')
  print(f'-- Running:')
  print(f'--   {cmddis.pop()}')

  while len(cmddis):      #  evenutally pretty-ify this to put options and args on same line,
    opt  = cmddis.pop()   #  but need care to figure out stuff like '-threads 4 /input/file'
    optl = len(opt) + 1   #  vs '-inputs /input/files'

    print(f'--     {opt}')

  #print(cmdarg)
  #print(' '.join(cmdarg))

  print(f'--  stdin  < {stdin}')
  print(f'--  stdout > {stdout}')
  print(f'--  stderr > {stderr}')
  print(f'--')

  try:
    sp = subprocess.Popen(cmdarg, stdin=stdin, stdout=stdout, stderr=stderr)
  except OSError:
    print(f'ERROR: failed to run {command}: Not found.')
    #  Failed to find binary
    pass
  except ValueError:
    #  Invalid arguments to Popen
    print(f'ERROR: failed to run {command}: Invalid arguments.')
    pass
  except subprocess.TimeoutExpired:
    #  Ran too long; not expected since we're not asking for a timeout.
    pass

  return sp

