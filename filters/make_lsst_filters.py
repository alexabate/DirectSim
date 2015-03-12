"""
Basic script to take LSST baseline throughputs and write them to DirectSim
in the correct format

usage: make_lsst_filters.py [location_of_baseline_filters]

e.g.

make_lsst_filters.py ~/Software/throughputs/baseline/

"""
import sys

thruputs_loc = sys.argv[1]
filters = ["u","g","r","i","z","y"]


for filt in filters:

    print "On filter",filt

    # open filter file
    f = open(thruputs_loc + "/total_" + filt + ".dat", 'r')
    
    # first read header line: assumes there is only one
    line =  f.readline()

    # new file to write to
    fnew = open(filt + "_lsst_throughputs.txt", 'w')

    # then read rest of lines
    #ii = 0
    while line:
        line = f.readline()
        #if ii<3:
        #    print float(line.split(' ')[0])
        #ii+=1
        if len(line)>0:
            wl_meters = float(line.split(' ')[0])*1e-9
            trans = float(line.split(' ')[1])
            fnew.write(str(wl_meters) + "  " + str(trans) + "\n")
 
    # then close both files    
    fnew.close()
    f.close()


