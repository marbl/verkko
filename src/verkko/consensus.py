import sys
import inspect

#  See comments in verkko.py, please.
def desc():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''Compute contig consensus sequences'''
  return inspect.cleandoc(s)

def help():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''Consensus ...

        Nulla non risus ac magna venenatis scelerisque. Nam et aliquet
        justo. Nam accumsan dui eu risus aliquam, ut varius tortor
        auctor. Pellentesque bibendum, ligula ac sollicitudin viverra, lacus
        elit elementum nulla, quis venenatis sem justo efficitur
        ligula. Nullam rutrum nisi non odio aliquam, congue vestibulum magna
        placerat. Vestibulum ante ipsum primis in faucibus orci luctus et
        ultrices posuere cubilia curae; Fusce eleifend vehicula
        faucibus. Proin non lectus sit amet justo efficitur dictum. Maecenas
        ullamcorper mattis metus vitae posuere. Aenean sollicitudin lectus
        vitae dui accumsan, eu scelerisque nulla accumsan. Nam posuere magna
        vel augue euismod, vel aliquet augue vestibulum. Quisque tincidunt
        est ante, non venenatis nunc egestas fermentum. Maecenas facilisis
        lorem sit amet velit bibendum, eget tempus lorem vehicula. Aenean et
        ante est.'''
  return inspect.cleandoc(s)


def computeConsensus(args):
  print(args)

