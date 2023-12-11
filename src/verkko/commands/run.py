import sys
import inspect

import importlib

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



#  See comments in verkko.py, please.
def synopsis():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''run all steps of Verkko on the local machine'''
  return inspect.cleandoc(s)

def details():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''run ...

        Vestibulum ut gravida metus, nec laoreet massa. Integer vitae
        ullamcorper felis, pellentesque aliquam nisi. Interdum et malesuada
        fames ac ante ipsum primis in faucibus. Donec ac dui sed mi vehicula
        aliquam non non metus. Nulla mattis augue sapien, eget aliquet urna
        tincidunt vitae. Phasellus dignissim est sollicitudin nunc dignissim
        accumsan vel id est. Phasellus vel arcu placerat, convallis neque eu,
        vestibulum mi. Morbi blandit et enim eget hendrerit. Nulla
        facilisi. Quisque ornare elementum quam, nec bibendum elit molestie
        id. Proin in enim eu dolor accumsan dapibus non in odio. Nulla ut
        mauris ut felis ullamcorper accumsan sed non turpis. Cras viverra sit
        amet velit in consequat. Cras efficitur arcu nisi, id ornare odio
        viverra ut. In hac habitasse platea dictumst.'''
  return inspect.cleandoc(s)

def run(args):
  loadHiFi.run(args)
  loadONT.run(args)
  loadParental.run(args)
  loadHiC.run(args)

  correctHiFi.run(args)
  correctONT.run(args)
  correctParental.run(args)
  correctHiC.run(args)

  buildGraph.run(args)
  refineGraph.run(args)

  buildPaths.run(args)
  refinePaths.run(args)

  haplotype.run(args)
  phase.run(args)
  rukki.run(args)

  consensus.run(args)
