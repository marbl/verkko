import os
import sys
import inspect

import argparse
import configparser

import importlib

theory          = importlib.import_module('commands.theory')

run             = importlib.import_module('commands.run')          #  Slightly wonky import method so we can use
dryrun          = importlib.import_module('commands.dry-run')      #  module files named after the subcommand.

loadHiFi        = importlib.import_module('commands.load-hifi-reads')
loadONT         = importlib.import_module('commands.load-ont-reads')
loadParental    = importlib.import_module('commands.load-parental-reads')
loadHiC         = importlib.import_module('commands.load-hic-reads')

correctHiFi     = importlib.import_module('commands.correct-hifi-reads')
correctONT      = importlib.import_module('commands.correct-ont-reads')
correctParental = importlib.import_module('commands.correct-parental-reads')
correctHiC      = importlib.import_module('commands.correct-hic-reads')

buildGraph      = importlib.import_module('commands.build-hifi-graph')
refineGraph     = importlib.import_module('commands.refine-graph')

buildPaths      = importlib.import_module('commands.build-ont-paths')
refinePaths     = importlib.import_module('commands.refine-paths')

haplotype       = importlib.import_module('commands.haplotype')
phase           = importlib.import_module('commands.phase')
rukki           = importlib.import_module('commands.rukki')

consensus       = importlib.import_module('commands.generate-consensus')


#  Description and main help text.
#    shortdesc is displayed when verkko is run with no arguments, and for '-h' and '--help'.
#    longhelp  is displayed when for the 'help' subcommand.
#
#  For modules:
#    synopsis() - A (very) short description, printed when run with no args or '-h', etc.
#    details()  - Long pre-formatted text for 'subcommand -h'
#
def shortdesc():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''A hybrid genome assembly pipeline for telomere-to-telomere assembly of PacBio HiFi
        and Oxford Nanopore reads.'''
  return inspect.cleandoc(s)




#  Create a master argument parser and suppress displaying the -h option.
#  The lambda ets the formatter to a reasonable max width.
#
#  https://docs.python.org/2/library/argparse.html#argparse.RawTextHelpFormatter
#  We use the RawTextHelpFormatter so we can pre-format description strings.
#   -- but what is the width=77 parameter??
#    HelpFormatter               -
#    RawTextHelpFormatter        -
#    RawDescriptionHelpFormatter - do not reformat 'description' or 'epilog'

pFmt = lambda prog: argparse.RawTextHelpFormatter(prog, width=77)
pMain = argparse.ArgumentParser(prog=sys.argv[0], add_help=False, description=shortdesc(), formatter_class=pFmt, epilog='https://github.com/marbl/verkko')
pMain.add_argument('-h', '--help', action='help', help=argparse.SUPPRESS)

#  Add another parser to use for holding options that are shared among all
#  the sub-parsers, then add a (suppressed) help option and any global
#  options.

pParent = argparse.ArgumentParser(add_help=False)
pParent.add_argument('-h', '--help', action='help', help=argparse.SUPPRESS)

pParent.add_argument('-d', '--directory',      dest='workdir',      default=None,    required=True,                       help='Directory where Verkko will compute the assembly.')

#Parent.add_argument('--keep-intermediate',    dest='saveinter',    default=True,    action='store_const', const=True,    help='')
pParent.add_argument('--discard-intermediate', dest='saveinter',    default=True,    action='store_const', const=False,   help='remove intermediate files as soon as possible')
pParent.add_argument('--cleanup',              dest='cleanup',      default=False,   action='store_const', const=True,    help='remove all intermediate files when assembly is complete')
#Parent.add_argument('--no-cleanup',           dest='cleanup',      default=False,   action='store_const', const=False,   help='')

pParent.add_argument('--python',               dest='python',       default=None,                                         help='path to a python interpreter')
pParent.add_argument('--perl',                 dest='perl',         default=None,                                         help='path to a perl interpreter')

pParent.add_argument('--mbg',                  dest='mbg',          default=None,                                         help='path to the MBG executable')
pParent.add_argument('--graphaligner',         dest='graphaligner', default=None,                                         help='path to the GraphAligner executable')
pParent.add_argument('--mashmap',              dest='mashmap',      default=None,                                         help='path to the mashmap executable')
pParent.add_argument('--winnowmap',            dest='winnowmap',    default=None,                                         help='path to the WinnowMap executable')
pParent.add_argument('--bwa',                  dest='bwa',          default=None,                                         help='path to the BWA executable')
pParent.add_argument('--samtools',             dest='samtools',     default=None,                                         help='path to the samtools exeutable')

pParent.add_argument('--local',                dest='grid',         default='local', action='store_const', const='local', help='run compute jobs on the local machine')
pParent.add_argument('--sge',                  dest='grid',         default='local', action='store_const', const='sge',   help='run compute jobs using SGE')
pParent.add_argument('--slurm',                dest='grid',         default='local', action='store_const', const='slurm', help='run compute jobs using Slurm')
pParent.add_argument('--lsf',                  dest='grid',         default='local', action='store_const', const='lsf',   help='run compute jobs using LSF')

pParent.add_argument('--local-memory',         dest='local_mem',    help='set memory limit for --local')
pParent.add_argument('--local-cpus',           dest='local_cpus',   help='set CPU limit for --local')

pParent.add_argument('--snakeopts',            dest='', help='')

#
#  Create a sub-command parser for all the sub-commands.  Do NOT make this
#  'required=True', as that will suppress the long help message when we're
#  run with no sub-command.
#
#  Generate sub-parsers for each command.
#    help        is displayed in the global help message.  Shoud be short.
#    description is displayed during in sub-command help.  Can be long.
#

pCmds = pMain.add_subparsers(title='Commands', dest='subcommand')
pC    = {}

pC['theory']                 = pCmds.add_parser('theory',              help=theory.synopsis(),          description=theory.details(),          add_help=False, formatter_class=pFmt)

pC['run']                    = pCmds.add_parser('run',                 help=run.synopsis(),             description=run.details(),             add_help=False, formatter_class=pFmt, parents=[pParent])
pC['dry-run']                = pCmds.add_parser('dry-run',             help=dryrun.synopsis(),          description=dryrun.details(),          add_help=False, formatter_class=pFmt, parents=[pParent])

pC['load-hifi-reads']        = pCmds.add_parser('load-hifi-reads',     help=loadHiFi.synopsis(),        description=loadHiFi.details(),        add_help=False, formatter_class=pFmt, parents=[pParent])
pC['load-ont-reads']         = pCmds.add_parser('load-ont-reads',      help=loadONT.synopsis(),         description=loadONT.details(),         add_help=False, formatter_class=pFmt, parents=[pParent])
pC['load-parental-reads']    = pCmds.add_parser('load-parental-reads', help=loadParental.synopsis(),    description=loadParental.details(),    add_help=False, formatter_class=pFmt, parents=[pParent])
pC['load-hic-reads']         = pCmds.add_parser('load-hic-reads',      help=loadHiC.synopsis(),         description=loadHiC.details(),         add_help=False, formatter_class=pFmt, parents=[pParent])

pC['correct-hifi-reads']     = pCmds.add_parser('correct-hifi-reads',  help=correctHiFi.synopsis(),     description=correctHiFi.details(),     add_help=False, formatter_class=pFmt, parents=[pParent])
pC['correct-ont-reads']      = pCmds.add_parser('correct-ont-reads',   help=correctONT.synopsis(),      description=correctONT.details(),      add_help=False, formatter_class=pFmt, parents=[pParent])
pC['correct-trio-reads']     = pCmds.add_parser('correct-trio-reads',  help=correctParental.synopsis(), description=correctParental.details(), add_help=False, formatter_class=pFmt, parents=[pParent])
pC['correct-hic-reads']      = pCmds.add_parser('correct-hic-reads',   help=correctHiC.synopsis(),      description=correctHiC.details(),      add_help=False, formatter_class=pFmt, parents=[pParent])

pC['build-hifi-graph']       = pCmds.add_parser('build-hifi-graph',    help=buildGraph.synopsis(),      description=buildGraph.details(),      add_help=False, formatter_class=pFmt, parents=[pParent])
pC['refine-graph']           = pCmds.add_parser('refine-hifi-graph',   help=refineGraph.synopsis(),     description=refineGraph.details(),     add_help=False, formatter_class=pFmt, parents=[pParent])

pC['build-ont-paths']        = pCmds.add_parser('align-ont-to-graph',  help=buildPaths.synopsis(),      description=buildPaths.details(),      add_help=False, formatter_class=pFmt, parents=[pParent])
pC['refine-paths']           = pCmds.add_parser('refine-paths',        help=refinePaths.synopsis(),     description=refinePaths.details(),     add_help=False, formatter_class=pFmt, parents=[pParent])

pC['compute-haplotype']      = pCmds.add_parser('compute-haplotype',   help=haplotype.synopsis(),       description=haplotype.details(),       add_help=False, formatter_class=pFmt, parents=[pParent])
pC['compute-phase']          = pCmds.add_parser('compute-phase',       help=phase.synopsis(),           description=phase.details(),           add_help=False, formatter_class=pFmt, parents=[pParent])
pC['rukki']                  = pCmds.add_parser('rukki',               help=rukki.synopsis(),           description=rukki.details(),           add_help=False, formatter_class=pFmt, parents=[pParent])

pC['generate-consensus']     = pCmds.add_parser('generate-consensus',  help=consensus.synopsis(),       description=consensus.details(),       add_help=False, formatter_class=pFmt, parents=[pParent])

#
#  Add options and defaults to each sub-command
#

pC['theory'].set_defaults(func=theory.run)

pC['run'].set_defaults(func=run.run)
pC['dry-run'].set_defaults(func=dryrun.run)

pC['load-hifi-reads'].set_defaults(func=loadHiFi.run)
pC['load-ont-reads'].set_defaults(func=loadONT.run)
pC['load-parental-reads'].set_defaults(func=loadParental.run)
pC['load-hic-reads'].set_defaults(func=loadHiC.run)

pC['correct-hifi-reads'].set_defaults(func=correctHiFi.run)
pC['correct-ont-reads'].set_defaults(func=correctONT.run)
pC['correct-trio-reads'].set_defaults(func=correctParental.run)
pC['correct-hic-reads'].set_defaults(func=correctHiC.run)

pC['build-hifi-graph'].set_defaults(func=buildGraph.run)
pC['refine-graph'].set_defaults(func=refineGraph.run)

pC['build-ont-paths'].set_defaults(func=buildPaths.run)
pC['refine-paths'].set_defaults(func=refinePaths.run)

pC['compute-haplotype'].set_defaults(func=haplotype.run)
pC['compute-phase'].set_defaults(func=phase.run)
pC['rukki'].set_defaults(func=rukki.run)

pC['generate-consensus'].set_defaults(func=consensus.run)

#
#  Load defaults.
#    From the system  verkko.ini (as found via the path)
#    From the system  verkko.ini (as found via the environment)
#    From the user    verkko.ini (as found via the environment)
#    From the current verkko.ini (in the current directory)
#    From the command line
#

def loadConfig(config, dirname, filename):
  if dirname != None:
    filename = os.path.join(dirname, filename)
  if filename == None:
    return
  if not os.path.exists(filename):
    return

  #print(f"Loading config from '{filename}'.")

  cIn = configparser.ConfigParser(strict=False)
  cIn.read(filename)

  for s in cIn:
    for k in cIn[s]:
      if s not in cConfig:
        config[s] = {}
      config[s][k] = cIn[s][k]

cConfig   = configparser.ConfigParser(strict=False)

rundir = os.path.dirname(os.path.abspath(__file__))
sysdir = os.getenv('VERKKO_HOME')
usrdir = os.getenv('HOME')
cwddir = os.getcwd()

loadConfig(cConfig, rundir, 'verkko.ini')
loadConfig(cConfig, sysdir, 'verkko.ini')
loadConfig(cConfig, usrdir, 'verkko.ini')
loadConfig(cConfig, cwddir, 'verkko.ini')

#cConfig.write(sys.stdout)

#  Parse the args, show help, or run the command.
#  If there is no subcommand, show global help, listing the subcommands.
#  If the subcommand is 'help', show the theory text above.
#  If there is no workdir, show help for the subcommand.

pArgs = pMain.parse_args(sys.argv[1:])

if pArgs.subcommand == None:
  pMain.print_help()

#elif pArgs.subcommand == 'help':
#  pC[pArgs.subcommand].print_help()

#elif pArgs.workdir == None:
#  pC[pArgs.subcommand].print_help()

else:
  pArgs.func(pMain)

exit(0)
