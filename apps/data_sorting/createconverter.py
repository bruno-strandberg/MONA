#!/usr/bin/python
"""
This script auto-generates a reader class for a ECAP PID summary file in ROOT format and auto-generates a code skeleton for an application that uses the auto-generated reader to convert data to MONA format. For example, executing `./createconverter.py -s my_input_file.root -n NewReader -t PID -i iRODS_location -c "application info"` will do the following:

1. Look for a `TTree` named "PID" in the file "my_input_file.root" and use the function `TTree::MakeClass("parsers/NewReader")` to auto-create a ROOT class "parsers/NewReader.h/C" for reading the PID tree.
2. Create a skeleton for an application `NewReader_to_MONA.C` that includes the reader class and sets up translation of the data to analysis format. The user will need to modify the application to decide which variables are mapped to `SummaryEvent`.
3. Use the iRODS location and (-i) and comment (-c) to automatically document the new application in README.md

After modification of the application `NewReader_to_MONA.C`, the user should type `make` and the application should be built.

NB! Because of some (Py)ROOT issue, the script says 
" Error in <TClass::LoadClassInfo>: no interpreter information for class THbookTree is available even though it has a TClass initialization routine. "
As long as the readers are created, this error message can be ignored.

Usage:
    createconverer.py -s SUMMARYFILE -n READERNAME [-t TREENAME] [-i IRODSLOCATION] [-c COMMENT]
    createconverer.py -h                                                                     
                                                                                          
Option:
    -s SUMMARYFILE              PID summary file from ECAP in .root format
    -n READERNAME               Name of the auto-generated reader, recommended `ECAPYYMMDD`
    -t TREENAME                 Name of the ROOT tree in SUMMARYFILE [default: PID]
    -i IRODSLOCATION            iRODS location of the PID output file the new reader and application operate on [default: ]
    -c COMMENT                  Comment that is used to automatically document the new application in `README.md`
    -h --help                   Show this screen

"""

from docopt import docopt
args = docopt(__doc__)
import os
from ROOT import *

#======================================================================================
# first thing, generate a ROOT reader class for the new PID file
#======================================================================================

fin = TFile(args['-s'], "READ")

if not fin.IsOpen():
    raise Exception("ERROR in createconverter.py: cannot open file {}".format(args['-s']) )

tin = fin.Get(args['-t'])

if tin == None:
    raise Exception("ERROR in createconverter.py: cannot find TTree {} in file {}".format(args['-t'], args['-s']) )

tin.MakeClass( "{}".format( args['-n'] ) )
os.system("mv {}.* ./parsers".format(args['-n']) )

#======================================================================================
# now create a skeleton application for conversion
#======================================================================================

appname = "{}_to_MONA.C".format(args['-n'])

conversion = { "NEWREADERCLASS":args['-n'] , "NEWPIDINPUTFILE":args['-s'], "NEWTREENAME":args['-t']}


template = open('skeleton/TemplateApp.C','r')
app = open(appname, 'w')

for line in template.readlines():

    linein  = line
    lineout = linein

    for key, value in conversion.iteritems():

        if key in line:
            lineout = linein.replace(key, value)
            linein  = lineout # this helps if two keys on one line

    app.write(lineout)

app.close()

#======================================================================================
# add the application description to the README.md
#======================================================================================

appdata = "\nApplication `{0}` and parsers `parsers/{1}.h/C` - comment: {2} ; iRODS location of the ECAP PID output that the application operates on: {3}\n".format(appname, args['-n'], args['-c'], args['-i'])

readme = open('README.md', 'a')
readme.write( appdata )
readme.close()
