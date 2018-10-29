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
__date__ = '$Date: 2017-12-01 21:16:30 (Fri, 01 Dec 2017) $'[7:-21]
__version__ = '$Version: copy_results.py V1.0 $'[10:-2]


################################################################################
# FUNCTION: to copy results from source directory to target directory
################################################################################
def copy_results():
    if len(sys.argv) == 6:
        yyyy = int(sys.argv[1])
        doy = int(sys.argv[2])
        ndays = int(sys.argv[3])
        dir_source = sys.argv[4]
        dir_target = sys.argv[5]
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
        dir_target_yyyyddd = os.path.join(dir_target,cyyyyddd)
        if not os.path.exists(dir_target_yyyyddd):
            os.makedirs(dir_target_yyyyddd)

        for dirname_source in glob.iglob(os.path.join(dir_source_yyyyddd,'result_*')):
            if os.path.isdir(dirname_source):
                print dirname_source
                dirname = os.path.basename(dirname_source)
                dirname_target = os.path.join(dir_target_yyyyddd,dirname)
                print dirname_target
                shutil.copytree(dirname_source,dirname_target)

        doy += 1


################################################################################
# Main program
################################################################################
if __name__ == '__main__':
    copy_results()