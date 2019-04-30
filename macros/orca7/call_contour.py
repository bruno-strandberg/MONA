#!/usr/bin/python

th23values         = [0.46, 0.58]
InvertedOrdered    = [True, False]
IncludeSystematics = [True, False]
IncludePriors      = [True, False]

for th23 in th23values:
    for IO in InvertedOrdered:
        for IS in IncludeSystematics:
            for IP in IncludePriors:
            
                syscmd = "./contour -t {}".format(th23)

                # choose ordering
                if IO:
                    syscmd += " -i"

                # include systematics or not
                if IS:
                    syscmd += " -s"

                # if using systematics, use priors or not
                if IP:
                    syscmd += " -p"

                print syscmd
