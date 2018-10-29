function sod=transT_hms2sod(ih,im,sec)
%SYNTAX:
%     =================================
%     | sod=transT_hms2sod(ih,im,sec) |
%     =================================
%
%     hour, minute, second to seconds of day
%
%INPUTS:
%     ih: hour
%     im: minute
%    sec: second
%
%OUTPUT:
%    sod: seconds of day
%
%originally written by Feng Zhou on 2014/09/22 @ East China Normal Univ.
%
%email: zhouforme@163.com
%
%Dept. of Information Science Technology, East China Normal Univ., China

%%===========================BEGIN PROGRAM====================================%%
sod=ih*3600.0+im*60.0+sec;

return;
%%===========================END PROGRAM======================================%%