import sys
import textwrap
import inspect
import argparse
import configparser

import correct
import mbg
import consensus


#  Description and main help text.
#    shortdesc is displayed when verkko is run with no arguments, and for '-h' and '--help'.
#    longhelp  is displayed when for the 'help' subcommand.
#
#  For modules:
#    desc() - A (very) short description, printed when run with no args or '-h', etc.
#    help() - Long pre-formatted text for 'subcommand -h'
#
def shortdesc():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''A hybrid genome assembly pipeline for telomere-to-telomere assembly of PacBio HiFi
        and Oxford Nanopore reads.'''
  return inspect.cleandoc(s)

def longhelp():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''Verkko is a hybrid genome assembly pipeline developed for telomere-to-telomere
        assembly of PacBio HiFi and Oxford Nanopore reads. Verkko is Finnish for net, mesh
        and graph.

        Verkko uses Canu to correct remaining errors in the HiFi reads, builds a multiplex
        de Bruijn graph using MBG, aligns the Oxford Nanopore reads to the graph using
        GraphAligner, progressively resolves loops and tangles first with the HiFi reads
        then with the aligned Oxford Nanopore reads, and finally creates contig consensus
        sequences using Canu's consensus module.

        Verkko will create a subdirectory to compute all results in, the -d option.
          verkko -d output-directory [...].'''
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
pParent.add_argument('-d', '--directory', dest='workdir', help='Directory where Verkko will compute the assembly.')

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

pC['help']      = pCmds.add_parser('help',      help=shortdesc(),      description=longhelp(),       add_help=False, formatter_class=pFmt)
pC['correct']   = pCmds.add_parser('correct',   help=correct.desc(),   description=correct.help(),   add_help=False, formatter_class=pFmt, parents=[pParent])
pC['mbg']       = pCmds.add_parser('mbg',       help=mbg.desc(),       description=mbg.help(),       add_help=False, formatter_class=pFmt, parents=[pParent])
pC['consensus'] = pCmds.add_parser('consensus', help=consensus.desc(), description=consensus.help(), add_help=False, formatter_class=pFmt, parents=[pParent])

#  Set default action for each command.

pC['help']     .set_defaults(func=help)
pC['correct']  .set_defaults(func=correct.correctHiFi)
pC['mbg']      .set_defaults(func=mbg.buildGraph)
pC['consensus'].set_defaults(func=consensus.computeConsensus)

#  Parse the args, show help, or run the command.
#  If there is no subcommand, show global help, listing the subcommands.
#  If the subcommand is 'help', show the theory text above.
#  If there is no workdir, show help for the subcommand.

pArgs = pMain.parse_args(sys.argv[1:])

if pArgs.subcommand == None:
  pMain.print_help()

elif pArgs.subcommand == 'help':
  pC[pArgs.subcommand].print_help()

elif pArgs.workdir == None:
  pC[pArgs.subcommand].print_help()

else:
  pArgs.func(pMain)

exit(0)
