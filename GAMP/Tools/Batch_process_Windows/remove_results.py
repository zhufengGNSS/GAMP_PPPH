#!/usr/bin/env python
# coding:utf-8

################################################################################
# PROGRAM: 
################################################################################



################################################################################
# Import Python modules
import os,sys,pydoc
import shutil,glob


__author__ = 'Feng Zhou @ ECNU'
__date__ = '$Date: 2017-12-03 19:54:30 (Sun, 03 Dec 2017) $'[7:-21]
__version__ = '$Version: remove_results.py V1.0 $'[10:-2]


################################################################################
# FUNCTION: to remove results from source directory
################################################################################
def remove_results():
    if len(sys.argv) == 5:
        yyyy = int(sys.argv[1])
        doy = int(sys.argv[2])
        ndays = int(sys.argv[3])
        dir_source = sys.argv[4]
    else:
        return 0

    cyyyy = str(yyyy)
    for i in range(ndays):
        if doy < 10:
            cdoy = '00'+str(doy)
        elif doy >= 10 and doy < 100:
            cdoy = '0'+str(doy)
        else:
            cdoy = str(doy)

        cyyyyddd = cyyyy + cdoy
        dir_source_yyyyddd = os.path.join(dir_source,cyyyyddd)
        for dirname_source in glob.iglob(os.path.join(dir_source_yyyyddd,'result_*')):
        #for dirname_source in glob.iglob(os.path.join(dir_source_yyyyddd,'result_GE_sta_DF_noGIM_wum')):
            if os.path.isdir(dirname_source):
                print dirname_source
                shutil.rmtree(dirname_source)

        doy += 1


################################################################################
# Main program
################################################################################
if __name__ == '__main__':
    remove_results()