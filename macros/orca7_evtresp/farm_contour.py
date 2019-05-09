#!/usr/bin/python

import os

nohup              = False
th23values         = [0.46, 0.58]
InvertedOrdered    = [True, False]
IncludeSystematics = [True, False]
IncludePriors      = [True, False]

#==============================================================================
# create contour execution commands
#==============================================================================

os.system( "mkdir {}/tmp-contour".format(os.getcwd()) )
syscmds = []

for th23 in th23values:
    for IO in InvertedOrdered:
        for IS in IncludeSystematics:
            for IP in IncludePriors:

                syscmd  = "{}/./contour -t {}".format( os.getcwd(), th23 )
                outname = "{}/tmp-contour/contour_th23={}".format( os.getcwd(), th23 )

                # choose ordering
                if IO:
                    syscmd += " -i"
                    outname += "_IO"
                else:
                    outname += "_NO"

                # include systematics or not
                if IS:
                    syscmd += " -s"
                    outname += "_WithSys"
                else:
                    outname += "_NoSys"

                # if using systematics, use priors or not
                if IP:
                    syscmd += " -p"
                    outname += "_WithPriors"
                else:
                    outname += "_NoPriors"

                if IP and not IS:
                    continue
                else:
                    syscmd += " -o {}.root".format(outname)
                    syscmds.append( syscmd )

#==============================================================================
# create scripts to farm the commands
#==============================================================================

jobfiles = []

for i, cmd in enumerate(syscmds):

    jobfilename = "{}/tmp-contour/job_{}.sh".format(os.getcwd(), i)

    jobfile = open(jobfilename, 'w')

    # copy lines from bashrc
    bashrcf = open("/user/bstrand/.bashrc", 'r')
    copylines = False
    for line in bashrcf:

        if "source" in line:
            copylines = True

        if copylines:
            jobfile.write(line)

    bashrcf.close()
    jobfile.write( cmd )
    jobfile.close()
    jobfiles.append( jobfilename )

#==============================================================================
# execute the commands
#==============================================================================

for job in jobfiles:
    
    if nohup:
        farmcmd = "nohup bash {} &".format(job)
    else:
        farmcmd = "qsub -q short7 -o {0} -e {0} {1}".format(os.getcwd()+"/tmp-contour", job)

    os.system(farmcmd)
