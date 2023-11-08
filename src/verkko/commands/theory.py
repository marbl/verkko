import sys
import inspect

#  See comments in verkko.py, please.
def synopsis():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''show theory of operation'''
  return inspect.cleandoc(s)

def details():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''theory-of-operation

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
          verkko -d output-directory [].'''
  print(inspect.cleandoc(s))
