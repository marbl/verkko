import sys
import inspect

#  See comments in verkko.py, please.
def desc():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''Generate Canu-corrected HiFi reads'''
  return inspect.cleandoc(s)

def help():
  #     -90-columns-------------------------------------------------------------------------------
  s=f'''Correction ...

        Praesent vel elementum velit. Phasellus vel sapien vehicula, mollis
        nulla malesuada, varius ante. Aenean porttitor erat enim, nec
        accumsan lacus convallis vel. Sed tempus dapibus consequat. Aenean
        posuere et elit nec aliquet. Nulla eu ullamcorper arcu, eget tempus
        ligula. Nulla tincidunt commodo nulla, quis semper ligula porta
        a. Nunc in augue vel justo lobortis molestie. Integer iaculis risus
        ut metus mattis, nec lacinia velit condimentum.'''
  return inspect.cleandoc(s)


def correctHiFi(args):
  print(args)

