% Change the delay on line 10 to change how much delay you would like to
% approximate.

clear all;
rng(1);

% Single joint reaching movements:
G = .14;        % Viscous Constant: Ns/m
I = .1;         % Inertia Kgm2
tau = 0.066;    % Muscle time constant, s
delay = .25; 

A1 = [0 1 0
    0 -G/I 1/I
    0 0 -1/tau];
B1 = [0 0 1/tau]';
C1 = [1 0 0];

A2 = [0 1 0 0
    0 -G/I 4/(delay*I) -1/I
    0 0 -2/delay 1
    0 0 0 -1/tau];
B2 = [0 0 0 1/tau]';
C2 = [1 0 0 0];

G = ss(A1,B1,C1,0);
H = ss(A2,B2,C2,0);
freq = 0.01:0.01:1000;
[mag1, phase1] = bode(G,freq);
[mag2, phase2] = bode(H,freq);
mag1 = 20*log10(mag1);
mag2 = 20*log10(mag2);
phase1 = squeeze(phase1);
phase2 = squeeze(phase2)-360;

lag = phase1 - phase2;
lag = 2*pi*1000*lag./(freq'*360);

figure
subplot(2,1,1)
semilogx(freq,phase1)
hold on
title('Open-Loop Plant')
semilogx(freq,phase2)
ylabel('Phase (deg)')
legend({'Undelayed',['Delay = ',num2str(delay*1000),' ms']})

subplot(2,1,2)
semilogx(freq,lag)
xlabel('Frequency (Hz)')
ylabel('Time Difference (ms)')
