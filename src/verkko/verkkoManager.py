import os
#rom   snakemake.io import glob_wildcards, expand
from   collections  import namedtuple, deque
import subprocess

import verkkoHelper  as vH    #  For 'ruleOutput' and 'ruleInputs'


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


class verkkoTask:
  def __init__(self, command, *args,
               stdin   = '/dev/null',
               stdout  = '/dev/null',
               stderr  = '/dev/null',
               threads = 1,
               memory  = 1024):

    if stdin  == None:   stdin  = '/dev/null'
    if stdout == None:   stdout = '/dev/null'
    if stderr == None:   stderr = '/dev/null'

    self.command = command    #  Path to program to run.
    self.args    = args       #  Arguments to the program.
    self.stdin   = stdin      #  -
    self.stdout  = stdout     #  - Input/Output filename or None/'/dev/null'.
    self.stderr  = stderr     #  -
    self.threads = threads    #  Number of CPUs that will be used.
    self.memory  = memory     #  Amount of memory, MB, that will be used.

    self.process = None       #  A subprocess.Popen() result.

  def display(self, label):
    lines  = []
    cmdarg = [self.command] + vH.flatten(self.args)
    cmddis = cmdarg.copy()
    cmddis.reverse()

    lines.append(f'----------------------------------------')
    lines.append(f'-- {label}:')
    lines.append(f'--   {cmddis.pop()}')

    while len(cmddis):      #  evenutally pretty-ify this to put options and args on same line,
      opt  = cmddis.pop()   #  but need care to figure out stuff like '-threads 4 /input/file'
      optl = len(opt) + 1   #  vs '-inputs /input/files'

      lines.append(f'--     {opt}')

    lines.append(f'--   < {self.stdin}')
    lines.append(f'--  1> {self.stdout}')
    lines.append(f'--  2> {self.stderr}')
    lines.append(f'--')

    return lines

  def execute(self):
    fstdin  = subprocess.DEVNULL
    fstdout = subprocess.DEVNULL
    fstderr = subprocess.DEVNULL

    if self.stdin  != '/dev/null':   fstdin  = open(self.stdin,  'r')
    if self.stdout != '/dev/null':   fstdout = open(self.stdout, 'w')
    if self.stderr != '/dev/null':   fstderr = open(self.stderr, 'w')

    cmdarg = [self.command] + vH.flatten(self.args)
    cmddis = cmdarg.copy()
    cmddis.reverse()

    for l in self.display('Running:'):
      print(l)

    try:
      sp = subprocess.Popen(cmdarg, stdin=fstdin, stdout=fstdout, stderr=fstderr)
    except OSError as e:
      raise vH.verkkoException(f'ERROR: failed to run {e.filename}: {e.strerror}.')
    except ValueError:
      raise vH.verkkoException(f'ERROR: failed to run {command}: Invalid arguments.')

    return sp





class verkkoManager:
  def __init__(self):
    self.processes = []

  def submit(self, command, *args, stdin=None, stdout=None, stderr=None, threads=1, memory=16):
    self.processes.append(verkkoTask(command, args,
                                     stdin=stdin, stdout=stdout, stderr=stderr,
                                     threads=threads, memory=memory))

  def execute(self):
    for c in self.processes:
      c.process = c.execute()

  def appendlog(self, logname):
    log = []

    if os.path.isfile(logname):
      lll = deque([], 50)

      log.append(f'--  {logname}:')

      try:
        with open(logname) as f:
          lll.extend(f)
      except:
        pass

      if len(lll) > 0:
        for l in lll:
          log.append('--    ' + l.rstrip())
      else:
        log.append(f'--    <empty>')

      log.append(f'--')

    return log


  def wait(self):
    failed = []

    for c in self.processes:
      if c.process:
        print(f'Waiting on \'{c.command}\'')
        r = c.process.wait()

        if (r != 0):
          failed.append('--------------------')

          if (r < 0):
            failed.append(f'Failed with signal {-r}:')
          if (r > 0):
            failed.append(f'Failed with return code {r}:')

          failed.extend(c.display('Failed:'))
          failed.extend(self.appendlog(c.stdout))
          failed.extend(self.appendlog(c.stderr))


    if len(failed):
      raise vH.verkkoException(failed)
