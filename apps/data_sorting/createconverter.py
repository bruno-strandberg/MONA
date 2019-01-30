#!/usr/bin/python
"""
This script auto-generates a reader class for a ECAP PID summary file in ROOT format and auto-generates a code skeleton for an application that uses the auto-generated reader to convert data to NNMO format. For example, executing `./createconverter.py -s my_input_file.root -n NewReader -t PID` will do the following:

1. Look for a `TTree` named "PID" in the file "my_input_file.root" and use the function `TTree::MakeClass("NewReader")` to auto-create a ROOT class "NewReader.h/C" for reading the PID tree.
2. Create a skeleton for an application `NewReaderToSummary.C` that includes the reader class and sets up translation of the data to analysis format. The user will need to modify the application to decide which variables are mapped to `SummaryEvent`.
3. Modify `Makefile` by adding `NewReader.C` to the list of common source files; after modifying the application `NewReaderToSummary.C`, the user should type `make` and the application should be built.

NB! Because of some (Py)ROOT issue, the script says 
" Error in <TClass::LoadClassInfo>: no interpreter information for class THbookTree is available even though it has a TClass initialization routine. "
As long as the readers are created, this error message can be ignored.

Usage:
    createconverer.py -s SUMMARYFILE -n READERNAME [-t TREENAME]
    createconverer.py -h                                                                     
                                                                                          
Option:
    -s SUMMARYFILE              PID summary file from ECAP in .root format
    -n READERNAME               Name of the auto-generated reader, e.g. PIDAlpha, PIDBeta, ...
    -t TREENAME                 Name of the ROOT tree in SUMMARYFILE [default: PID]
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

#======================================================================================
# now create a skeleton application for conversion
#======================================================================================

appname = "{}ToSummary.C".format(args['-n'])

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
# modify the makefile
#======================================================================================

make_out = open("Makefile.tmp",'w')
make_in  = open("Makefile",'r')

for line in make_in.readlines():

    lineout = line
    if "COMMON" in line and args['-n'] not in line and "$(COMMON)" not in line:
        lineout = line[ :line.rfind('C')+1 ] + " {}.C\n".format(args['-n'])

    make_out.write(lineout)

make_in.close()
make_out.close()
os.system( "mv {} {}".format("Makefile.tmp","Makefile") )
