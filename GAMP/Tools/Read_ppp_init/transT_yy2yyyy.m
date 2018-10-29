function yyyy=transT_yy2yyyy(yy)
%SYNTAX:
%     ===========================
%     | yyyy=transT_yy2yyyy(yy) |
%     ===========================
%
%     get 4-digit year, valid from 1951-2050
%
%INPUTS:
%     yy: 2- or 4-digit year
%
%OUTPUT:
%   yyyy: 4-digit year
%
%originally written by Feng Zhou on 2014/09/22 @ East China Normal Univ.
%
%email: zhouforme@163.com
%
%Dept. of Information Science Technology, East China Normal Univ., China

%%===========================BEGIN PROGRAM====================================%%
if (yy>1900)
    yyyy=yy;
    return;
elseif (yy>50&&yy<1900)
    yyyy=yy+1900;
elseif (yy<=50)
    yyyy=yy+2000;
end

return;
%%===========================END PROGRAM======================================%%