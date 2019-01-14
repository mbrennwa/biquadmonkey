% m-file program to determine biquad filter coefficients for an IIR filter from a given filter transfer function:
%
% - load target transfer function data from an ASCII data file (column-1: frequency in Hz, column-2: gain in dB, column-3: phase in degrees)
% - ask for FDLS parameter range (maximum numbers of poles, zeros, and time delay)
% - ask for DSP sampling rate
% - fit IIR filter to target function using FDLS method
% - determine biquad coefficients for IIR filter
% - plot transfer functions of IIR filter vs target filter
% - print biquad coefficients for use with miniDSP on console
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

% load data file with target response curve here (magnitude in dB and phase in degrees)
disp ('Select data file with target filter function (f/mag/phase)...')
[FNAME, FPATH, FLTIDX] = uigetfile ();
u = load ([FPATH FNAME]);

% ask for n and fs:
fs = abs (input('DSP sampling rate (Hz, default: fs=96000): '));
if isempty(fs)
	fs = 96000
end

NP  = round (input('Maximum number of poles allowed in the filter (default: NP=16): '));
if isempty(NP)
	NP = 16
end
if NP < 1
	NP = 1
end

NZ  = round (input('Maximum number of zeros allowed in the filter (default: NZ=16): '));
if isempty(NZ)
	NZ = 16
end
if NZ < 1
	NZ = 1;
end

NT = abs (input('Maximum time delay (number of samples) allowed to improve stability of the IIR filter (default: NT = 10): '));
if isempty(NT)
	NT = 10
end

disp('');
disp ('Fitting biquad filter to target...');

% Prepare data:
f_target = u(:,1);
mag_target = u(:,2);
phase_target = u(:,3);

fLow  = min (f_target);
fHigh = max (f_target);

% remove target data above Nyquist frequency of DSP sampling rate:
k = find (f_target < fs/2);
f_target = f_target(k);
mag_target = mag_target(k);
phase_target = phase_target(k);

% convert target response (magnitude/dB and phase/deg) to complex transfer function H_target(z):
w   = pi*f_target/(fs/2); % normalized frequency (0...pi)
r   = 10.^(mag_target / 20); % magnitude of H_target
phi = phase_target / 180*pi; % phase of H_target
H_target = r .* exp(i*phi);

% fit IIR filter to the target transfer function using all possible parameter combinations

% wgt = 1./sqrt(abs(H_target));
wgt = 1./abs(H_target).^0.4;
x = [];
k = 1;
for np = 1:NP
	for nz = 1:NZ for nt = 1:NT
		[B,A] = fdls (H_target,w,np,nz,fs,wgt,nt); % determine B and A for current np,nz,dt
		H = freqz (B,A,w); % evaluate transfer function of current B and A
		dH = H - H_target; % residuals to target
		d = sum ( abs(dH).^2 ./ wgt ); % weighted sum of residuals
		x = [ x ; [ np nz nt d ] ];
		k = k+1;
	end end
	disp (sprintf('Progress: %g %%',np/NP*100))
end

% find best fit:
[u,k]=min(x(:,4));
np = x(k,1);
nz = x(k,2);
nt = x(k,3);
disp ('');
disp (sprintf('Best fit result for NP = %i poles, NZ = %i zeros and NT = %i delay samples.',np,nz,nt));
disp ('');

% re-calculate B and A:
[B,A] = fdls (H_target,w,np,nz,fs,wgt,nt);

% plot poles and zeros in the z-plane
figure(1);
clf;
zplane(B(:)',A(:)');
title ('Poles (x) and zeros (0)')

% evaluate transfer function of the fitted IIR filter
f_IIR         = logspace(log10(fLow/2),log10(fs/2),1000);
[H_IIR,W_IIR] = freqz (B,A,2*pi*f_IIR/fs);
mag_IIR       = 20*log10(abs(H_IIR));
phase_IIR     = arg(H_IIR)/pi*180;
f_IIR         =  W_IIR/pi/2*fs;

% determine poles and zeros of IIR filter:
[Z,P,G] = tf2zp (B,A);

% check for stability:
if any (abs(P) >= 1)
	% Not all poles are inside the unit circle!
	disp (''); disp ('');
	disp ('***** WARNING: the IIR filter is not stable!!! *****')
end

% convert poles and zeros to biquads:
BIQ = zp2sos (Z,P,G);

% plot results
figure(2);
clf;
k_target = find ( f_target > 0 ); k_IIR = find ( f_IIR > 0 );
ff = [ min([f_target(k_target)(:);f_IIR(k_IIR)(:)]) fs/2 ];
subplot(2,1,1);
semilogx ( f_target(k_target),mag_target(k_target),'r-','linewidth',2 , f_IIR(k_IIR),mag_IIR(k_IIR),'k-','linewidth',1 );
grid on
xlim(ff);
ylabel ('Gain (dB)')
legend ('Target','Biquad')
subplot(2,1,2);
semilogx ( f_target(k_target),phase_target(k_target),'r-','linewidth',2 , f_IIR(k_IIR),phase_IIR(k_IIR),'k-','linewidth',1 );
grid on
xlim(ff);
ylabel ('Phase (deg.)')
xlabel ('Frequency (Hz)');


% print miniDSP coefficients to console:
n = max([np,nz]);
disp('');
disp ('***** miniDSP biquad coefficients')
disp (sprintf('***** sampling rate: %g Hz',fs))
disp (sprintf('***** filter order: %i',n))
disp (sprintf('***** number of zeros: %i',nz))
disp (sprintf('***** number of poles: %i',np))
disp (sprintf('***** time delay: %i samples',nt))

disp('');

N_BIQ = N_BIQ = max ( [ 8 , round(n/2) ] );
if N_BIQ > 8
	disp (sprintf('***** WARNING: filter order n=%i. This requires %i biquads, but the miniDSP only allows 8 per DSP !!!',n,N_BIQ))
end

for k = 1:N_BIQ
	if k <= size(BIQ,1)
		a0 =  1;
		a1 = -BIQ(k,5) / BIQ(k,4);
		a2 = -BIQ(k,6) / BIQ(k,4);
		b0 =  BIQ(k,1) / BIQ(k,4);
		b1 =  BIQ(k,2) / BIQ(k,4);
		b2 =  BIQ(k,3) / BIQ(k,4);

	else
		a0 = 1;
		a1 = 0;
		a2 = 0;
		b0 = 1;
		b1 = 0;
		b2 = 0;
		
	end

	disp (sprintf('biquad%i,',k))
	disp (sprintf('b0=%.15f,',b0))
	disp (sprintf('b1=%.15f,',b1))
	disp (sprintf('b2=%.15f,',b2))
	disp (sprintf('a1=%.15f,',a1))
	if k < N_BIQ
		disp (sprintf('a2=%.15f,',a2))
	else
		disp (sprintf('a2=%.15f',a2))
	end
end
