#!/usr/bin/env python
# coding:utf-8


################################################################################
# PROGRAM: 
################################################################################
"""

 Batch processing in precise point positioning (PPP) mode using GAMP software

 Usage: python gamp_batch.py yyyy doy ndays ac satsys mode freq ion

 OPTIONS:
   <yyyy>         - 4-digit year
   <doy>          - Day of year
   <ndays>        - Number of consecutive days of data to retrieve
   <ac>           - GNSS Analysis Center, it can be igs/cod/esa/gfz/grg/com/gbm/grm/wum
   <satsys>       - Satellite system, e.g., G, R, C, GR, GC, or GRCE
   <mode>         - Processing mode, e.g., sta for static PPP, and kin for kinematic PPP
   <freq>         - Frequency selection, i.e., SF for single-frequency, and DF for dual-frequency
   <ion>          - Y denotes using ionospheric-constrain, N denotes not using

EXAMPLES: python gamp_batch.py 2017 001 1 com G kin DF Y
          python gamp_batch.py 2017 001 1 grm GR kin SF N
          python gamp_batch.py 2017 001 1 gbm GRCE sta SF N
          python gamp_batch.py 2017 001 1 wum GRCE kin DF N

Changes: 22-Nov-2017   fzhou: create the prototype of the scripts for version 1.0

 to get help, type:
           python gamp_batch.py
 or        python gamp_batch.py -h

"""
################################################################################
# Import Python modules
import os,sys,pydoc,types
import shutil,platform
import numpy as np
from daytime import *
from utility import *


__author__ = 'Feng Zhou @ ECNU'
__date__ = '$Date: 2017-11-22 16:52:30 (Wed, 22 Nov 2017) $'[7:-21]
__version__ = '$Version: gamp_batch.py V1.0 $'[10:-2]


def gamp_batch():
    if '-h' in sys.argv or len(sys.argv) < 2:
        pydoc.help(os.path.basename(sys.argv[0]).split('.')[0])
        return 0
    elif (len(sys.argv) == 9):
        yyyy = int(sys.argv[1])
        doy = int(sys.argv[2])
        ndays = int(sys.argv[3])
        ac = sys.argv[4]
        satsys = sys.argv[5]
        mode = sys.argv[6]
        freq = sys.argv[7]
        ion = sys.argv[8]
    else:
        return 0

    if 'Windows' in platform.system():
        dir_path = 'C:\mannual_GAMP\Examples'
        dir_data = 'C:\mannual_GAMP\Examples'
    elif 'Linux' or 'Mac' in platform.system():
        dir_path = '/home/fzhou/data'
        dir_data = '/home/fzhou/data'

    dir_table = os.path.join(dir_path,'tables')
    dir_igs = os.path.join(dir_data,'igs_prod')
    dir_snx = os.path.join(dir_data,'igs_snx')
    dir_dcb = os.path.join(dir_data,'dcb')
    dir_ion = os.path.join(dir_data,'ion')
    dir_nav = os.path.join(dir_data,'nav')
    cfgfil = os.path.join(dir_table,'gamp.cfg')
    ocnfil = os.path.join(dir_table,'ocnload.blq')
    crdfil = os.path.join(dir_table,'site.crd')

    cyyyy = str(yyyy)
    cyy = cyyyy[2:4]
    for i in range(ndays):
        if doy < 10:
            cdoy = '00'+str(doy)
        elif doy >= 10 and doy < 100:
            cdoy = '0'+str(doy)
        else:
            cdoy = str(doy)

        cyyyyddd = cyyyy + cdoy
        datdir_local = os.path.join(dir_path,cyyyyddd)
        if not os.path.exists(datdir_local):
            os.makedirs(datdir_local)

        dir_res = os.path.join(datdir_local,'result')
        if os.path.exists(dir_res):
            shutil.rmtree(dir_res)

        doy1 = doy - 1
        doy2 = doy + 1
        gpsw, weekd = yrdoy2gpst(yyyy, doy)
        gpsw1, weekd1 = yrdoy2gpst(yyyy, doy1)
        gpsw2, weekd2 = yrdoy2gpst(yyyy, doy2)
        mm, dd = yrdoy2ymd(yyyy, doy)
        if mm < 10:
            cmm = '0'+str(mm)
        elif mm >= 10 and mm <= 12:
            cmm = str(mm)

        # replace configure file and others
        filename = os.path.join(datdir_local,'gamp.cfg')
        if os.path.isfile(filename):
            os.remove(filename)
        shutil.copy(cfgfil,filename)
        filename = os.path.join(datdir_local,'ocnload.blq')
        if not os.path.isfile(filename):
            shutil.copy(ocnfil,filename)
        filename = os.path.join(datdir_local,'site.crd')
        if not os.path.isfile(filename):
            shutil.copy(crdfil,filename)
        if gpsw >= 1934:
            filename = os.path.join(datdir_local,'igs14.atx')
            if os.path.isfile(filename):
                os.remove(filename)
            if ac == 'com':
                atxfil = os.path.join(dir_table,'igs14_igs.atx')
                shutil.copy(atxfil,filename)
            elif ac == 'gbm':
                atxfil = os.path.join(dir_table,'igs14_gbm.atx')
                shutil.copy(atxfil,filename)
            elif ac == 'grm':
                atxfil = os.path.join(dir_table,'igs14_igs.atx')
                shutil.copy(atxfil,filename)
            elif ac == 'wum':
                atxfil = os.path.join(dir_table,'igs14_wum.atx')
                shutil.copy(atxfil,filename)
        else:
            filename = os.path.join(datdir_local,'igs08.atx')
            if os.path.isfile(filename):
                os.remove(filename)
            if ac == 'com':
                atxfil = os.path.join(dir_table,'igs08_igs.atx')
                shutil.copy(atxfil,filename)
            elif ac == 'gbm':
                atxfil = os.path.join(dir_table,'igs08_gbm.atx')
                shutil.copy(atxfil,filename)
            elif ac == 'grm':
                atxfil = os.path.join(dir_table,'igs08_igs.atx')
                shutil.copy(atxfil,filename)
            elif ac == 'wum':
                atxfil = os.path.join(dir_table,'igs08_wum.atx')
                shutil.copy(atxfil,filename)

        # replace sp3 and clk files
        ext = '*.sp3'
        removefils(datdir_local,ext)
        ext = '*.clk'
        removefils(datdir_local,ext)
        
        prefix = [ac + str(gpsw1) + str(weekd1), ac + str(gpsw) + str(weekd), ac + str(gpsw2) + str(weekd2)]
        suffix = ['.sp3', '.clk']
        for pre in prefix:
            for suf in suffix:
                fil = pre + suf
                fil_source = os.path.join(dir_igs,fil)
                fil_target = os.path.join(datdir_local,fil)
                shutil.copy(fil_source,fil_target)

        #copy sinex file
        snxf = 'igs' + str(gpsw) + '.snx'
        snx_source = os.path.join(dir_snx,snxf)
        snx_target = os.path.join(datdir_local,snxf)
        if not os.path.isfile(snx_target):
            if not os.path.isfile(snx_source):
                print '\n ***ERROR: SNX file NOT found ',snx_source
                return
            else:
                shutil.copy(snx_source,snx_target)

        #copy eop file
        eopf = 'igs' + str(gpsw) + '7.erp'
        eop_source = os.path.join(dir_igs,eopf)
        eop_target = os.path.join(datdir_local,eopf)
        if not os.path.isfile(eop_target):
            if not os.path.isfile(eop_source):
                print '\n ***ERROR: EOP file NOT found ',eop_source
                return
            else:
                shutil.copy(eop_source,eop_target)

        #copy dcb files
        p1p2_dcb = 'P1P2' + cyy + cmm + '.DCB'
        p1c1_dcb = 'P1C1' + cyy + cmm + '.DCB'
        p2c2_dcb = 'P2C2' + cyy + cmm + '.DCB'
        mgex_dcb = 'CAS0MGXRAP_' + cyyyy + cdoy + '0000_01D_01D_DCB.BSX'
        dcbfils = [p1p2_dcb,p1c1_dcb,p2c2_dcb,mgex_dcb]
        for dcbf in dcbfils:
            dcb_source = os.path.join(dir_dcb,dcbf)
            dcb_target = os.path.join(datdir_local,dcbf)
            if not os.path.isfile(dcb_target):
                if not os.path.isfile(dcb_source):
                    print '\n ***ERROR: DCB file NOT found ',dcb_source
                    return
                else:
                    shutil.copy(dcb_source,dcb_target)

        #copy ion file
        ionf = 'CODG' + cdoy + '0.' + cyy + 'I'
        ion_source = os.path.join(dir_ion,ionf)
        ion_target = os.path.join(datdir_local,ionf)
        if not os.path.isfile(ion_target):
            if not os.path.isfile(ion_source):
                print '\n ***ERROR: ION file NOT found ',ion_source
                return
            else:
                shutil.copy(ion_source,ion_target)

        #copy nav file
        navf = 'brdm' + cdoy + '0.' + cyy + 'p'
        nav_source = os.path.join(dir_nav,navf)
        nav_target = os.path.join(datdir_local,navf)
        if not os.path.isfile(nav_target):
            if not os.path.isfile(nav_source):
                print '\n ***ERROR: NAV file NOT found ',nav_source
                return
            else:
                shutil.copy(nav_source,nav_target)

        ext = '*.' + cyy + 'o'
        nobs = numfiles(datdir_local,ext)
        if nobs == 0:
            dir_mgex = os.path.join(dir_data,'mgex_obs')
            dir_obs = os.path.join(dir_mgex,cyyyyddd)
            site_list = os.path.join(dir_path,'site.list')
            sitl = np.loadtxt(site_list, dtype=np.str)
            sitlst = sitl.tolist()
            if type(sitlst) == types.ListType:
                for site in sitlst:
                    site = site.lower()
                    ofil = site + cdoy + '0.' + cyy + 'o'
                    obs_source = os.path.join(dir_obs,ofil)
                    obs_target = os.path.join(datdir_local,ofil)
                    shutil.copy(obs_source,obs_target)
            elif type(sitlst) == types.StringType:
                site = sitlst.lower()
                ofil = site + cdoy + '0.' + cyy + 'o'
                obs_source = os.path.join(dir_obs,ofil)
                obs_target = os.path.join(datdir_local,ofil)
                shutil.copy(obs_source,obs_target)

        # modified the configure file
        filename = os.path.join(datdir_local,'gamp.cfg')
        ssss = '                    = ' + datdir_local + '\n'
        string_switch(filename,'C:\mannual_GAMP\Examples',ssss)
        ssss = 'start_time          = 0    2017/01/01  00:00:00.0       %(0:from obs  1:from inp)\n'
        string_switch(filename,'start_time',ssss)
        ssss = 'end_time            = 0    2017/01/01  00:00:00.0       %(0:from obs  1:from inp)\n'
        string_switch(filename,'end_time',ssss)
        if satsys == 'G':
            ssss = 'navsys              = 1                     %(1:gps  4:glo  5:gps+glo  8:gal  16:qzs  32:bds)\n'
            string_switch(filename,'navsys',ssss)
        elif satsys == 'R':
            ssss = 'navsys              = 4                     %(1:gps  4:glo  5:gps+glo  8:gal  16:qzs  32:bds)\n'
            string_switch(filename,'navsys',ssss)
            ssss = 'gloicb              = 0                     %(0:off  1:linear function of frequency number  2:quadratic polynomial function of frequency number\n'
            string_switch(filename,'gloicb',ssss)
        elif satsys == 'C':
            ssss = 'navsys              = 32                    %(1:gps  4:glo  5:gps+glo  8:gal  16:qzs  32:bds)\n'
            string_switch(filename,'navsys',ssss)
        elif satsys == 'GR':
            ssss = 'navsys              = 5                     %(1:gps  4:glo  5:gps+glo  8:gal  16:qzs  32:bds)\n'
            string_switch(filename,'navsys',ssss)
            ssss = 'gloicb              = 0                     %(0:off  1:linear function of frequency number  2:quadratic polynomial function of frequency number\n'
            string_switch(filename,'gloicb',ssss)
        elif satsys == 'GC':
            ssss = 'navsys              = 33                    %(1:gps  4:glo  5:gps+glo  8:gal  16:qzs  32:bds)\n'
            string_switch(filename,'navsys',ssss)
        elif satsys == 'GE':
            ssss = 'navsys              = 9                     %(1:gps  4:glo  5:gps+glo  8:gal  16:qzs  32:bds)\n'
            string_switch(filename,'navsys',ssss)
        elif satsys == 'GJ':
            ssss = 'navsys              = 17                    %(1:gps  4:glo  5:gps+glo  8:gal  16:qzs  32:bds)\n'
            string_switch(filename,'navsys',ssss)
        elif satsys == 'GRC':
            ssss = 'navsys              = 37                    %(1:gps  4:glo  5:gps+glo  8:gal  16:qzs  32:bds)\n'
            string_switch(filename,'navsys',ssss)
            ssss = 'gloicb              = 0                     %(0:off  1:linear function of frequency number  2:quadratic polynomial function of frequency number\n'
            string_switch(filename,'gloicb',ssss)
        elif satsys == 'GRCE':
            ssss = 'navsys              = 45                    %(1:gps  4:glo  5:gps+glo  8:gal  16:qzs  32:bds)\n'
            string_switch(filename,'navsys',ssss)
            ssss = 'gloicb              = 0                     %(0:off  1:linear function of frequency number  2:quadratic polynomial function of frequency number\n'
            string_switch(filename,'gloicb',ssss)
        elif satsys == 'GRCEJ':
            ssss = 'navsys              = 61                    %(1:gps  4:glo  5:gps+glo  8:gal  16:qzs  32:bds)\n'
            string_switch(filename,'navsys',ssss)
            ssss = 'gloicb              = 0                     %(0:off  1:linear function of frequency number  2:quadratic polynomial function of frequency number\n'
            string_switch(filename,'gloicb',ssss)

        if mode == 'sta':
            ssss = 'posmode             = 7                     %(0:spp  6:ppp_kinematic  7:ppp_static)\n'
            string_switch(filename,'posmode',ssss)
        elif mode == 'kin':
            ssss = 'posmode             = 6                     %(0:spp  6:ppp_kinematic  7:ppp_static)\n'
            string_switch(filename,'posmode',ssss)

        if freq == 'DF':
            ssss = 'inpfrq              = 2                     % the number of processing frequencies\n'
            string_switch(filename,'inpfrq',ssss)
            ssss = 'ionoopt             = 4                     %(0:off  1:brdc  2:IF12  3:UC1  4:UC12  5:ion-tec)\n'
            string_switch(filename,'ionoopt',ssss)
            ssss = 'prcNoise(ION)       = 4.0e-02               %process noise of random walk for slant ionospheric delay parameters (m/sqrt(s))\n'
            string_switch(filename,'prcNoise(ION)',ssss)
        elif freq == 'SF':
            ssss = 'inpfrq              = 1                     % the number of processing frequencies\n'
            string_switch(filename,'inpfrq',ssss)
            ssss = 'ionoopt             = 3                     %(0:off  1:brdc  2:IF12  3:UC1  4:UC12  5:ion-tec)\n'
            string_switch(filename,'ionoopt',ssss)
            ssss = 'prcNoise(ION)       = 4.0e-02               %process noise of random walk for slant ionospheric delay parameters (m/sqrt(s))\n'
            string_switch(filename,'prcNoise(ION)',ssss)

        if ion == 'Y':
            ssss = 'ionconstraint       = 1                     %(0:off  1:on)\n'
            string_switch(filename,'ionconstraint ',ssss)
        elif ion == 'N':
            ssss = 'ionconstraint       = 0                     %(0:off  1:on)\n'
            string_switch(filename,'ionconstraint ',ssss)

        cmd = 'gamp_c ' + filename
        os.system(cmd)

        dir_res = os.path.join(datdir_local,'result')
        if ion == 'Y':
            res_new = 'result_' + satsys + '_' + mode + '_' + freq + '_GIM_' + ac
            dir_res_new = os.path.join(datdir_local,res_new)
            os.rename(dir_res,dir_res_new)
        elif ion == 'N':
            res_new = 'result_' + satsys + '_' + mode + '_' + freq + '_noGIM_' + ac
            dir_res_new = os.path.join(datdir_local,res_new)
            os.rename(dir_res,dir_res_new)

        doy += 1


################################################################################
# Main program
################################################################################
if __name__ == '__main__':
    gamp_batch()