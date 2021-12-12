#############################################################################
 #
 #  This file is part of Verkko, a software program that assembles
 #  whole-genome sequencing reads into telomere-to-telomere
 #  haplotype-resolved chromosomes.
 #
 #  Except as indicated otherwise, this is a 'United States Government
 #  Work', and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
 #
 ##

import glob


##########
#
#  Paths to Verkko componenets.
#
#    VERKKO - the root directory of your Verkko installation.
#
#    PYTHON - name of yourpreffered python interpreter, or empty to use
#             whatever is in your environment.
#
VERKKO = config.get('VERKKO', 'please-set-VERKKO-in-config.yml')
PYTHON = config.get('PYTHON', '')


##########
#
#  Inputs!
#
#  While the shell should expand vars/globs in the input reads, silly
#  advanced users might add those in a hand-generated config file, and
#  so we need to expand wildcards ourself with:
#
#    HIFI_READS  = glob.glob(os.path.expandvars(config['HIFI_READS']))
#    ONT_READS   = glob.glob(os.path.expandvars(config['ONT_READS']))
#
#  But from verkko.sh, we get lists of expanded filenames directly.
#
HIFI_READS = config.get('HIFI_READS')
ONT_READS  = config.get('ONT_READS')


##########
#
#  Outputs!
#
#  The primary output is 'consensus'.
#
rule verkko:
    input:
        graph     = '5-untip/unitig-popped-unitig-normal-connected-tip.gfa',
        layout    = '6-layoutContigs/unitig-popped.layout',
        consensus = '7-consensus/unitig-popped.fasta'


##########
#
#  Local rules.
#
#  If your head node allows I/O tasks, include some of the following tasks
#  as local, to run them on the head node directly, instead of submitting
#  a job to the grid:
#
#    splitONT         - partitions ONT reads for alignment to the
#                       initial graph.
#
#    combineONT       - collects results of ONT-to-graph alignments.
#
#    buildPackages    - collects reads for contig consensus.  uses
#                       lots of memory but does no compute, except for
#                       sequence compression/decompression.
#
#    combineConsensus - collects results of contig consensus computations.
#
localrules: verkko


##########
#
#  Rules.
#
include: 'Snakefiles/functions.sm'

include: 'Snakefiles/1-buildGraph.sm'

include: 'Snakefiles/2-processGraph.sm'

include: 'Snakefiles/3-splitONT.sm'
include: 'Snakefiles/3-alignONT.sm'
include: 'Snakefiles/3-combineONT.sm'

include: 'Snakefiles/4-processONT.sm'

include: 'Snakefiles/5-untip.sm'

include: 'Snakefiles/6-layoutContigs.sm'

include: 'Snakefiles/7-extractONT.sm'
include: 'Snakefiles/7-buildPackage.sm'
include: 'Snakefiles/7-generateConsensus.sm'
include: 'Snakefiles/7-combineConsensus.sm'
