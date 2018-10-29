function y = nanrms(x,varargin)
%NANRMS	Root mean square.
%
% For vectors, NANRMS(x) returns the root mean square.
% For matrices, NANRMS(X) is a row vector containing the
% root mean square of each column. All NaN valuez
% are neglected.
%
% NANRMS(X,DIM) takes the rms along the dimension DIM of X. 
%
%See also: RMS, nanMEAN, nanMAX, nanMIN, STD

% $Id$
% $Date$
% $Author$
% $Revision$
% $HeadURL$
% $Keywords$

if isempty(x)
	y = NaN;
	return
end

    y = sqrt(nanmean(x.^2,varargin{:}));
