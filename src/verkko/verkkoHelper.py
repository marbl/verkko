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

  print(cmdarg)

  exit(0)

  try:
    sp = subprocess.Popen(cmdarg, stdin=stdin, stdout=stdout, stderr=stderr)
  except OSError:
    #  Failed to find binary
    pass
  except ValueError:
    #  Invalid arguments to Popen
    pass
  except subprocess.TimeoutExpired:
    #  Ran too long
    pass

  return sp

