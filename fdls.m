function [B,A] = fdls (H,w,nP,nZ,fs,wgt,nT)

% function [B,A] = fdls (H,w,nP,nZ,fs,wgt,nT)
% function fdls ('demo')
%
% Fit a digital filter to the target transfer function H(w) using the frequency-domain least-squares (FDLS) method. The digital filter transfer function is of the form H(z) = B(z)/A(z), where B and A are polynomials of degree nZ and nP.
% Additional weighting factors may be used with the fit. Also, a time delay may be used to improve the stability of the resulting filter.
%
% This code was heavily inspired by the FDLS_201.m file by Greg Berchin, who developed the FDLS method.
%
% Full descrition of the FDLS method: "Precise Filter Design" by Greg Berchin, doi: 10.1002/9780470170090.ch7, chapter 7 of "Streamlining Digital Signal Processing: A Tricks of the Trade Guidebook").
%
% [B,A] = fdls (H,w,nP,nZ,fs,wgt,nT) fits the digital filter with the INPUT and OUTPUT arguments as described below.
% 
% fdls ('demo') fits the digital filter as described in the section "Numerical Example" in the book chapter as cited above, and prints the X, Y, A and B coefficients. Note that the values in the 6th row of X and Y were computed wrong in the numerical example given in the book! The fdls.m function has been tested to give the same results as the original FDLS_201.m program by Greg Berchin.
%
% INPUT:
% H: complex transfer function
% w: frequency, normalised to sampling rate fs (radians, w = 0...pi)
% nP: number of poles of the digital filter
% nZ: number of zeros of the difital filter
% fs: sampling rate of the DSP (Hz)
% wgt (optional): weighting factors (dimensionless). Default: wgt = repmat (1,size(H))
% nT (optional): additional time delay of signal (number of samples). Default: nT = 0.
%
% OUTPUT (digital filter coefficients):
% B: coefficients of the nominator polynomial
% A: coefficients of the denominator polynomial
%
% DISCLAIMER:
% This file is part of BIQUADMONKEY.
% 
% BIQUADMONKEY is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% 
% BIQUADMONKEY is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with MATAA; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
% 
% Copyright (C) 2019 Matthias S. Brennwald.
% Contact: mbrennwa@gmail.com


% check if "demo" mode:
if nargin == 1
	if strcmp(H,'demo')
		isdemo = true;
		fM   = [ 0 19.685 35.4331 51.1811 59.0551 66.9291 106.299 389.764 ];
		magM = [ 0.2172 0.2065 0.1696 0.0164 1.3959 0.6734 0.349 0.3095 ];
		phiM = [ 0 -0.0156 -0.0383 3.0125 2.3087 0.955 0.0343 0.0031 ];
		fs  = 1000;
		H = magM .* exp(i*phiM);
		w = 2*pi * fM/fs;
		nP = 2; nZ = 2;
	else
		error('fdls: cannot parse input arguments.')
	end
else
	isdemo = false;
end


% make sure input arguments are right:
H = H(:);
w = w(:);
if any(w > pi)
	error ('fdls: data with frequency higher than Nyquist frequency is not allowed!')
end
if ~exist('wgt','var')
	wgt = [];
end
if isempty(wgt)
	wgt = repmat(1,size(w));
end
if ~exist('nT','var')
	nT = 0;
end
if nP < 0
	error ('fdls: number of poles must not be negative!')
end
if nZ < 0
	error ('fdls: number of zeros must not be negative!')
end

% convert / parse things:
f     = w/pi*(fs/2);	% natural frequency values (Hz)
mag   = abs(H(:));	% mangitude of target function
phase = arg(H(:));	% phase of target function
M     = length (f);	% number of "measurements" in the target function

% Check if fit is overspecified:
if M < (nP+nZ+1)
	error('fdls: system is underspecified (too few target data points / too many unknowns).')
end

% Add artificial delay to phase
phase = phase - nT*w;

% Create X matrix:
X = zeros(M,nP+nZ+1);
for m = 1:M
	X(m,nP+[1:nZ+1]) = cos(-[0:nZ] * w(m)) * wgt(m);
	if nP > 0 % number of poles > 0
		X(m,1:nP) = -mag(m) * cos((-[1:nP] * w(m)) + phase(m)) * wgt(m);
	end
end

% Create Y vector:
Y = mag .* cos(phase) .* wgt;

% Compute Theta:
Theta = X \ Y;

% B and A coefficients:
B = Theta((nP+1):(nP+1+nZ));
if (nP > 0)
    A = [ 1 ; Theta(1:nP) ];
else
    A = 1;
end

if isdemo
	disp(''); disp ('fm values:'); disp (fM(:));
	disp(''); disp ('Am values:'); disp (magM(:));
	disp(''); disp ('phim values:'); disp (phiM(:));
	disp(''); disp ('Y vector:'); disp (Y);
	disp(''); disp ('X matrix:'); disp (X);
	disp(''); disp ('B coefficients:'); disp (B);
	disp(''); disp ('A coefficients:'); disp (A);
	A = B = [];
end
