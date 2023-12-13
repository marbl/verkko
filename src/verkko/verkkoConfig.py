import os
import sys
import re
import argparse

#
#  Parse a standard .ini file, allowing:
#   - preservation of comments
#   - layering of input files
#   - command line override
#   - sections consisting entirely of filenames
#
#  Comment preservation works by keeping a copy of the input
#  lines in the order they are read.
#
#  Layering works by replacing the value of a section.parameter
#  as it is read from a subsequent file.
#
#  Command line overrides work by providing three subclasses of
#  argparse.Action that replace argparse actions 'store',
#  'store_const' and 'append'/'extend'.
#
#  File name sections are an extension to the ini format that
#  contain a list of files, one per line, instead of the usual
#  key=value format.
#



#
#  The intent of verkkoConfigLine was originally to parse out the pieces of
#  each line, then reconstruct the line with any updated values on output
#  while keeping the comment formatting the same.  This became infeasible
#  when parsing multiple ini files was allowed, because (a) comments could be
#  incorrect by later changes and (b) multi-line comments at the end of a
#  line would get intermixed.
#
#  --file 1--
#      [A]       #
#      val1 = 3  #  Be sure to set val1 and val2 to
#      val2 = 3  #  the same value, else Bad!!
#
#  --file 2--
#      [A]
#      val2 = 6    #  use more val2
#  
#  --merged result--
#      [A]       #
#      val1 = 3  #  Set val1 and val2 to
#      val2 = 6    #  use more val2
#
#  In the end, none of the fancy stuff was implemented.
#
class verkkoConfigLine(object):
  def __init__(self, indent, param, lgap, eq, section, rgap, value, cgap, comment):
    self.indent    = indent  #  {indent} {param}{lgap}{eq}{rgap}{value} {cgap}{#comment}
    self.param     = param   #  {indent} [{lgap}{section}{rgap}]        {cgap}{#comment}
    self.lgap      = lgap    #  {indent} {value}                        {cgap}{#comment}
    self.eq        = eq      #  {indent}                                      {#comment}
    self.section   = section
    self.rgap      = rgap
    self.value     = value
    self.valueorig = value
    self.cgap      = cgap
    self.comment   = comment if comment != None else ''

  def __str__(self):
    if   self.section and not self.param:
      return f'{self.indent}[{self.lgap}{self.section}{self.rgap}]{self.cgap}{self.comment}'
    elif self.param:
      return f'{self.indent}{self.param}{self.lgap}{self.eq}{self.rgap}{self.value}{self.cgap}{self.comment}'
    elif self.value:
      return f'{self.indent}{self.value}{self.cgap}{self.comment}'
    elif self.comment:
      return f'{self.indent}{self.comment}'
    else:
      return f''

#
#  Holds a map from string 'param-name' to string 'param-value'
#   - param-name MUST NOT be an integer.
#
#  Holds a list of name-less 'values' (e.g., input filenames)
#   - Access via integers.
#   - Access via assigning to a list
#
class verkkoConfigSection(object):
  def __init__(self, sectionname):
    self.sectionname = sectionname  #  String name of section
    self.params    = {}             #  Map from string param-name to string param-value
    self.values    = []             #  List of values with no param-name

  def __getitem__(self, param):     #  Returns string param-value for param-name.
    if   param in self.params:      #  If param-name is integer and not a param-name,
      return self.params[param]     #  returns the corresponding item.

    print(f'__getitem__ param {param} len {len(self.values)}')

    try:
      if int(param) < len(self.values):
        return self.values[int(param)]
    except:
      pass

    return None

  def __setitem__(self, param, value):
    self.params[param] = value




class verkkoConfig(object):
  def __init__(self):
    self.filenames = []   #  List of files we've loaded (to prevent loading same file again)
    self.sections  = {}   #  Map from section-name to 'map from param to value'
    self.lines     = []   #  List of verkkoConfigLine objects that created this verkkoConfig

    #  While all the configurable parameters/values come from the default
    #  verkko.ini file, it is VASTLY easier if we hardcode the input file
    #  sections and create them initially.
    #   * - easier in that we need to create these sections somewhere (else
    #       we blow up trying to iterate over the 'values' in the non-existent
    #       section) and so what better place to create than in the constructor?
    #
    for section in [ 'HIFI', 'ONT', 'HIC1', 'HIC2', 'HAP1', 'HAP2' ]:
      self.sections[section] = verkkoConfigSection(section)


  def __getitem__(self, section):     #  Return either the map of params -> values
    if section in self.sections:      #  for the specified section or an empty map.
      return self.sections[section]
    else:
      return {}

  def __iter__(self):
    return iter(self.sections)


  def load(self, dirname, filename):
    if filename == None:   return                  #  Silently ignore missing files.
    if dirname != None:                            #  Append dirname/filename.
      filename = os.path.join(dirname, filename)   #
    filename = os.path.abspath(filename)           #  Then find the full absolute path.

    if filename in self.filenames:    return       #  Don't load duplicates
    if not os.path.exists(filename):  return       #  Don't load missing files
    self.filenames.append(filename)                #  Remember.

    print(f"Loading config from '{filename}'.")

    self.lines.append(verkkoConfigLine('', None, None, None, None, None, None, None, f''))
    self.lines.append(verkkoConfigLine('', None, None, None, None, None, None, None, f'##########'))
    self.lines.append(verkkoConfigLine('', None, None, None, None, None, None, None, f'#'))
    self.lines.append(verkkoConfigLine('', None, None, None, None, None, None, None, f'#  Config from \'{filename}\''))
    self.lines.append(verkkoConfigLine('', None, None, None, None, None, None, None, f'#'))

    currentSection = ''

    with open(filename, 'rt') as f:
      for l in f:
        l = l.strip()

        comM = re.match(r'(\s*)(#.*?){0,1}\s*$', l)
        secM = re.match(r'(\s*)\[(\s*)(.+?)(\s*)](\s*)(#.*?){0,1}\s*$', l)
        parM = re.match(r'(\s*)(.+?)(\s*)([:=])(\s*)(.*?)(\s*)(#.*?){0,1}\s*$', l)
        valM = re.match(r'(\s*)(.+?)(\s*)(#.*?){0,1}\s*$', l)

        #  A comment-only line was parsed.
        if comM:
          indent  = comM[1]
          comment = comM[2]
          self.lines.append(verkkoConfigLine(indent, None, None, None, None, None, None, None, comment))

        #  A section header line was parsed.
        elif secM:
          indent  = secM[1]
          lgap    = secM[2]
          section = secM[3]
          rgap    = secM[4]
          cgap    = secM[5]
          comment = secM[6]
          self.lines.append(verkkoConfigLine(indent, None, lgap, None, section, rgap, None, cgap, comment))

          currentSection = section
          if currentSection not in self.sections:
            self.sections[currentSection] = verkkoConfigSection(section)

        #  A param=value line was parsed.
        elif parM:
          indent  = parM[1]
          param   = parM[2]
          lgap    = parM[3]
          eq      = parM[4]
          rgap    = parM[5]
          value   = parM[6]
          cgap    = parM[7]
          comment = parM[8]
          self.lines.append(verkkoConfigLine(indent, param, lgap, eq, currentSection, rgap, value, cgap, comment))

          #print(f'SET [{currentSection}][{param}] = \'{value}\'')
          self.sections[currentSection].params[param] = value

        #  A value-only line was parsed.
        elif valM:
          indent  = valM[1]
          value   = valM[2]
          cgap    = valM[3]
          comment = valM[4]
          self.lines.append(verkkoConfigLine(indent, None, None, None, None, None, value, cgap, comment))

          self.sections[currentSection].values.append(value)

        #  No idea what this line is supposed to be.
        else:
          print(f'Invalid line "{l}"')

        #print(f'LINE:  "{self.lines[-1]}"')

    self.lines.append(verkkoConfigLine('', None, None, None, None, None, None, None, f'#'))
    self.lines.append(verkkoConfigLine('', None, None, None, None, None, None, None, f'#  End of config from \'{filename}\''))
    self.lines.append(verkkoConfigLine('', None, None, None, None, None, None, None, f'#'))
    self.lines.append(verkkoConfigLine('', None, None, None, None, None, None, None, f'##########'))
    self.lines.append(verkkoConfigLine('', None, None, None, None, None, None, None, f''))


  def show(self):
    print(f'SHOW {self.lines}')
    for l in self.lines:
      print(l)

  def save(self, dirname, filename):
    if dirname != None:
      filename = os.path.join(dirname, filename)
    if filename == None:
      return

    print(f'Save config to \'{filename}\'')

    with open(filename, 'wt') as f:
      modified = False

      for l in self.lines:            #  Write all lines as they
        print(l, file=f)              #  appeared in the input.

      print(f'##########', file=f)
      print(f'#', file=f)
      print(f'#  Locally modified parameters', file=f)
      print(f'#', file=f)
      for l in self.lines:
        if l.param and l.value != self.sections[l.section][l.param]:
          l.value = self.sections[l.section][l.param]
          print(l, file=f)
      print(f'#', file=f)
      print(f'#  End of locally modified parameters', file=f)
      print(f'#', file=f)
      print(f'##########', file=f)
      print(f'', file=f)
      print(f'##########', file=f)
      print(f'#', file=f)
      print(f'#  Locally modified values', file=f)
      print(f'#', file=f)
      for sectionname in self.sections:
        sect = self.sections[sectionname]

        if len(sect.values) > 0:
          print(f'[{sectionname}]', file=f)
          for i in sect.values:
            print(i, file=f)
          print(f'', file=f)
      print(f'#', file=f)
      print(f'#  End of locally modified values', file=f)
      print(f'#', file=f)
      print(f'##########', file=f)



    print(f'Saved')




########################################
#
#  Override the argparse actions so we can store directly into the
#  proper section in a verkkoConfig.
#
#    derived from https://stackoverflow.com/a/16414640
#
class setopt(argparse.Action):
  def __init__(self, option_strings, vconfig, optsection, *args, **kwargs):
    self._vconfig    = vconfig
    self._optsection = optsection
    super(setopt, self).__init__(option_strings=option_strings, *args, **kwargs)

  def __call__(self, parser, namespace, value, option_string=None):
    self._vconfig[self._optsection][self.dest] = value

class setval(argparse.Action):
  def __init__(self, option_strings, vconfig, optsection, *args, **kwargs):
    self._vconfig    = vconfig
    self._optsection = optsection
    super(setval, self).__init__(option_strings=option_strings, *args, **kwargs, nargs=0)

  def __call__(self, parser, namespace, value, option_string=None):
    self._vconfig[self._optsection][self.dest] = self.const

class optapp(argparse.Action):
  def __init__(self, option_strings, vconfig, optsection, *args, **kwargs):
    self._vconfig    = vconfig
    self._optsection = optsection
    super(optapp, self).__init__(option_strings=option_strings, *args, **kwargs)

  def __call__(self, parser, namespace, values, option_string=None):
    #  All sections should exist already.
    #if self._optsection not in self._vconfig.sections:
    #  self._vconfig.sections[self._optsection] = verkkoConfigSection(self._optsection)
    self._vconfig[self._optsection].values += values



if __name__ == '__main__':
  vC = verkkoConfig()
  vC.load('.', 'verkko.ini')

  vC['MBG']['mbg-baseK'] = '1111'
  vC['HIC']['number-cpus'] = '888'

  vC.save('.', 'verkko.current.ini')
