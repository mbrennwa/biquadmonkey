% create a data file with a 3rd order Butterworth low-pass filter transfer function for testing biquadmonkey

fs = 44100;	% sampling rate
n  = 2;		% filter order
fc = 1000;      % filter cut-off frequency

% Butteworth filter:
[b,a] = butter(n,fc/(fs/2));

% Evaluate frequency response:
f = logspace (1,log10(fs/2),500);
H = freqz (b,a,f,fs);
m = 20*log10(abs(H));
p = arg(H) / pi *180;

f = f(1:end-1);
m = m(1:end-1);
p = p(1:end-1);

% Plot frequency response
subplot (2,1,1)
semilogx (f,m)
subplot (2,1,2)
semilogx (f,p)

% Save to test data file:
x = [ f(:) m(:)  p(:) ];
save -ascii butter_test.txt x

