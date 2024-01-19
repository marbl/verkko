#!/usr/bin/env python

import os
import sys
import inspect
import argparse
import traceback

import verkkoConfig as vC
import verkkoHelper as vH


import importlib

theory          = importlib.import_module('commands.theory')       #  Slightly wonky import method so we can use
                                                                   #  module files named after the subcommand,
run             = importlib.import_module('commands.run')          #  in particular, so that the file namess can
dryrun          = importlib.import_module('commands.dry-run')      #  have dashes in them.

setup           = importlib.import_module('commands.setup')

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








#
#  Load default values.
#    From the system  verkko.ini (as found via the path)
#    From the system  verkko.ini (as found via the environment)
#    From the user    verkko.ini (as found via the environment)
#    From the current verkko.ini (in the current directory)
#    From the command line
#
cErrors   = []
cConfig   = vC.verkkoConfig()

cConfig.load(os.path.dirname(os.path.abspath(__file__)), 'verkko.ini')  #  verkko defaults
cConfig.load(os.getenv('VERKKO_HOME'),                   'verkko.ini')  #  system defaults
cConfig.load(os.getenv('HOME'),                          'verkko.ini')  #  user defaults
cConfig.load(os.getcwd(),                                'verkko.ini')  #  run defaults

cConfig.lib = os.path.dirname(os.path.abspath(__file__)) + '/lib/verkko/'
cConfig.bin = os.path.dirname(os.path.abspath(__file__)) + '/lib/verkko/bin/'
cConfig.lib = '/work/verkko/lib/verkko/'
cConfig.bin = '/work/verkko/lib/verkko/bin/'

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

pParent.add_argument('--hifi',                                                          action=vC.optapp, nargs='+', vconfig=cConfig, optsection='HIFI')
pParent.add_argument('--ont', '--nano',                                                 action=vC.optapp, nargs='+', vconfig=cConfig, optsection='ONT')
pParent.add_argument('--hic1',                                                          action=vC.optapp, nargs='+', vconfig=cConfig, optsection='HIC1')
pParent.add_argument('--hic2',                                                          action=vC.optapp, nargs='+', vconfig=cConfig, optsection='HIC2')
#Parent.add_argument('--hap-kmers',                                                     action=vC.optapp, nargs='+', vconfig=cConfig, optsection='HAP1')
#Parent.add_argument('--hap-kmers',                                                     action=vC.optapp, nargs='+', vconfig=cConfig, optsection='HAP2')

#Parent.add_argument('--no-ont',                     dest='with-ont', default=True, action='store_const', const=False)

#Parent.add_argument('--paths',                      dest='hic1',   action='extend', nargs='+')
#Parent.add_argument('--assembly',                   dest='hic1',   action='extend', nargs='+')

#  Read Loading options
pParent.add_argument('--ont-split-bases',            dest='ont-split-bases',            action=vC.setopt, vconfig=cConfig, optsection='ONT_SPLIT')
pParent.add_argument('--ont-split-reads',            dest='ont-split-reads',            action=vC.setopt, vconfig=cConfig, optsection='ONT_SPLIT')
pParent.add_argument('--ont-split-min-length',       dest='ont-split-min-length',       action=vC.setopt, vconfig=cConfig, optsection='ONT_SPLIT')

pParent.add_argument('--hic-split-bases',            dest='hic-split-bases',            action=vC.setopt, vconfig=cConfig, optsection='HIC_SPLIT')
pParent.add_argument('--hic-split-reads',            dest='hic-split-reads',            action=vC.setopt, vconfig=cConfig, optsection='HIC_SPLIT')
pParent.add_argument('--hic-split-min-length',       dest='hic-split-min-length',       action=vC.setopt, vconfig=cConfig, optsection='HIC_SPLIT')

#  Canu Correction options
pParent.add_argument('--correct-kmer-size',          dest='kmer-size',                  action=vC.setopt, vconfig=cConfig, optsection='CANU_READ_CORRECTION')
pParent.add_argument('--correct-kmer-threshold',     dest='kmer-threshold',             action=vC.setopt, vconfig=cConfig, optsection='CANU_READ_CORRECTION')
pParent.add_argument('--correct-kmer-filter',        dest='kmer-filter',                action=vC.setopt, vconfig=cConfig, optsection='CANU_READ_CORRECTION')

pParent.add_argument('--no-correction',              dest='enabled',                    action=vC.setval, vconfig=cConfig, optsection='CANU_READ_CORRECTION', const='False')
pParent.add_argument('--correct-min-read-length',    dest='min-read-length',            action=vC.setopt, vconfig=cConfig, optsection='CANU_READ_CORRECTION')
pParent.add_argument('--correct-min-overlap-length', dest='min-overlap-length',         action=vC.setopt, vconfig=cConfig, optsection='CANU_READ_CORRECTION')
pParent.add_argument('--correct-min-kmers',          dest='min-kmers',                  action=vC.setopt, vconfig=cConfig, optsection='CANU_READ_CORRECTION')
pParent.add_argument('--correct-hash-bits',          dest='hash-bits',                  action=vC.setopt, vconfig=cConfig, optsection='CANU_READ_CORRECTION')

#  MBG Graph Building options
pParent.add_argument('--mbg-baseK',                  dest='mbg-baseK',                  action=vC.setopt, vconfig=cConfig, optsection='MBG')
pParent.add_argument('--mbg-maxK',                   dest='mbg-maxK',                   action=vC.setopt, vconfig=cConfig, optsection='MBG')
pParent.add_argument('--mbg-window',                 dest='mbg-window',                 action=vC.setopt, vconfig=cConfig, optsection='MBG')
pParent.add_argument('--mbg-max-resolution',         dest='mbg-max-resolution',         action=vC.setopt, vconfig=cConfig, optsection='MBG')
pParent.add_argument('--mbg-hifi-coverage',          dest='mbg-hifi-coverage',          action=vC.setopt, vconfig=cConfig, optsection='MBG')
pParent.add_argument('--mbg-unitig-abundance',       dest='mbg-unitig-abundance',       action=vC.setopt, vconfig=cConfig, optsection='MBG')

#  GraphAligner ONT Alignment options
pParent.add_argument('--ali-mxm-length',             dest='ali-mxm-length',             action=vC.setopt, vconfig=cConfig, optsection='ONT_ALIGN')
pParent.add_argument('--ali-mem-count',              dest='ali-mem-count',              action=vC.setopt, vconfig=cConfig, optsection='ONT_ALIGN')
pParent.add_argument('--ali-bandwidth',              dest='ali-bandwidth',              action=vC.setopt, vconfig=cConfig, optsection='ONT_ALIGN')
pParent.add_argument('--ali-multi-score-f',          dest='ali-multi-score-f',          action=vC.setopt, vconfig=cConfig, optsection='ONT_ALIGN')
pParent.add_argument('--ali-clipping',               dest='ali-clipping',               action=vC.setopt, vconfig=cConfig, optsection='ONT_ALIGN')
pParent.add_argument('--ali-min-score',              dest='ali-min-score',              action=vC.setopt, vconfig=cConfig, optsection='ONT_ALIGN')
pParent.add_argument('--ali-end-clipping',           dest='ali-end-clipping',           action=vC.setopt, vconfig=cConfig, optsection='ONT_ALIGN')
pParent.add_argument('--ali-incompat-cutoff',        dest='ali-incompat-cutoff',        action=vC.setopt, vconfig=cConfig, optsection='ONT_ALIGN')
pParent.add_argument('--ali-max-trace',              dest='ali-max-trace',              action=vC.setopt, vconfig=cConfig, optsection='ONT_ALIGN')
pParent.add_argument('--ali-seed-window',            dest='ali-seed-window',            action=vC.setopt, vconfig=cConfig, optsection='ONT_ALIGN')

#  Graph Refinement options
pParent.add_argument('--pop-min-allowed-cov',        dest='pop-min-allowed-cov',        action=vC.setopt, vconfig=cConfig, optsection='PROCESS_PATHS')
pParent.add_argument('--pop-resolve-steps',          dest='pop-resolve-steps',          action=vC.setopt, vconfig=cConfig, optsection='PROCESS_PATHS')
pParent.add_argument('--is-haploid',                 dest='is-haploid',                 action=vC.setval, vconfig=cConfig, optsection='PROCESS_PATHS', const='True')

#  Hi-C Phasing options
pParent.add_argument('--haplo-divergence',           dest='haplo-divergence',           action=vC.setopt, vconfig=cConfig, optsection='HIC')
pParent.add_argument('--uneven-depth',               dest='uneven-depth',               action=vC.setval, vconfig=cConfig, optsection='HIC', const='True')
pParent.add_argument('--no-rdna-tangle',             dest='no-rdna-tangle',             action=vC.setval, vconfig=cConfig, optsection='HIC', const='True')

#  Rukki options
pParent.add_argument('--ruk-enable',                 dest='ruk-enable',                 action=vC.setval, vconfig=cConfig, optsection='RUKKI', const='True')
pParent.add_argument('--ruk-hap1',                   dest='ruk-hap1',                   action=vC.setopt, vconfig=cConfig, optsection='RUKKI')
pParent.add_argument('--ruk-hap2',                   dest='ruk-hap2',                   action=vC.setopt, vconfig=cConfig, optsection='RUKKI')
pParent.add_argument('--ruk-type',                   dest='ruk-type',                   action=vC.setopt, vconfig=cConfig, optsection='RUKKI')
pParent.add_argument('--ruk-fract',                  dest='ruk-fract',                  action=vC.setopt, vconfig=cConfig, optsection='RUKKI')

#  Contig Filtering options
pParent.add_argument('--short-contig-length',        dest='short-contig-length',        action=vC.setopt, vconfig=cConfig, optsection='FILTERING')
pParent.add_argument('--screen',                     dest='screen',                     action=vC.setopt, vconfig=cConfig, optsection='FILTERING')



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

#C['setup']                  = pCmds.add_parser('setup',               help=setup.synopsis(),           description=setup.details(),           add_help=False, formatter_class=pFmt, parents=[pParent])

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


#  Parse the args, show help, or run the command.
#  If there is no subcommand, show global help, listing the subcommands.
#  If the subcommand is 'help', show the theory text above.
#  If there is no workdir, show help for the subcommand.

pArgs = pMain.parse_args(sys.argv[1:])

#  Copy a few choice parameters from pArgs to the config.


#  Dump the final config for later use.

cConfig.save(os.getcwd(), 'verkko.current.ini')

#  Sanitize the config.

setup.checkValid(cConfig, cErrors)      #  Are all parameters valid?
setup.checkInputs(cConfig, cErrors)     #  Do all input files exist?
setup.checkExternal(cConfig, cErrors)   #  Are all external tools available?
setup.checkConfig(cConfig, cErrors)     #  Are new config options compatible?
setup.failIfErrors(cConfig, cErrors)    #  Explode if any errors were detected.

#  Report usage if no sub-command.

#print(pArgs)


if pArgs.subcommand == None:
  pMain.print_help()

#elif pArgs.subcommand == 'help':
#  pC[pArgs.subcommand].print_help()

#elif pArgs.workdir == None:
#  pC[pArgs.subcommand].print_help()

else:
  try:
    vH.enterdir(pArgs.workdir)
    vH.leavedir()
    vH.enterdir(pArgs.workdir)
    vH.leavedir()
    vH.leavedir()

    pArgs.func(cConfig)

  except Exception as ve:
    print(f'')
    print(f'Verkko crashed.')
    print(f'')
    print(f'Traceback:')
    traceback.print_tb(ve.__traceback__)
    print(f'')
    print(f'Error:')
    print(f'  {ve}')
    print(f'')

exit(0)
