#!/usr/bin/env python
# coding:utf-8

################################################################################
# Import Python modules
import os
import glob


################################################################################
# FUNCTION: to find the specific files through a wildcard
################################################################################
def findfiles(rootdir,ext):
    result = []
    for filename in glob.iglob(os.path.join(rootdir,ext)):
        result.append(filename)
    return result


################################################################################
# FUNCTION: to count the specific files through a wildcard
################################################################################
def numfiles(rootdir,ext):
    nfil = 0
    for filename in glob.iglob(os.path.join(rootdir,ext)):
        nfil += 1
    return nfil


################################################################################
# FUNCTION: to remove the specific directories through a wildcard
################################################################################
def removedirs(rootdir,ext,show = False):
    for dirname in glob.iglob(os.path.join(rootdir,ext)):
        if show:
            print dirname
        os.removedirs(dirname)


################################################################################
# FUNCTION: to remove the specific files through a wildcard
################################################################################
def removefils(rootdir,ext,show = False):
    for filname in glob.iglob(os.path.join(rootdir,ext)):
        if show:
            print filname
        os.remove(filname)


################################################################################
# FUNCTION: to replace line-by-line in the file
################################################################################
def string_switch(filename,str1,str2,s=1):
    with open(filename, "r") as f:
        lines = f.readlines()
 
    with open(filename, "w") as f_w:
        n = 0
        if s == 1:
            for line in lines:
                if str1 in line:
                    line = line.replace(line,str2)
                    f_w.write(line)
                    n += 1
                    break
                f_w.write(line)
                n += 1
            for i in range(n,len(lines)):
                f_w.write(lines[i])
        elif s == 'g':
            for line in lines:
                if str1 in line:
                    line = line.replace(line,str2)
                f_w.write(line)