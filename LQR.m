% LQR model of sum-of-sines tracking. The first section of the code
% initializes all variables, simulates tracking, performs system ID. 
% All subsequent sections generate plots. The plots also compare simulation
% data to experimental data contained in 'dat.mat'. Free parameters to
% adjust are at the beginning of the first section (lag, Q, R).

clear all;
load('dat.mat');
rand('seed',1);

z = 1; % number of simulations
order = 3; % order of the system
delt = 0.005; % time step in secs

% adjust amount of time lag 
lag = 0.27/delt; % compute number of time steps of lag (0.27 = 270 milliseconds)

% values for Q and R taken from Qian infinite horizon model
Q = diag([1 0.01 0]); % accuracy cost- default is [1 0.01 0]
R = 0.0001; % effort cost- default is 0.0001

% parameters for A and B matrices
t1 = 0.224; 
t2 = 0.013;
t3 = 0.004;
k = 0;
b = t1 + t2;
m = t1*t2;
r = t3;

% generate A and B matrices in discrete time formulation
A = [0 1 0; -k/m -b/m 1/m; 0 0 -1/r];
A = eye(order) + delt*A;
B = [0 0 1/r]';
B = delt*B;

T = 42; % total simulation time

T2 = 40; % amount of analysis data
t = 0:delt:T2-delt; % x axis for graphs
nstep = round(T/delt); % number of simulation time steps
nstep2 = round(T2/delt); % number of analysis time steps

freqs_x = data.rot.avg.x_x.d.freqs; % frequencies of experimental x data
freqs_y = data.rot.avg.y_y.d.freqs; % frequencies of experimental y data

freq = (0.05:0.05:2.5)'; % frequencies used in the simulation
phases = 2*pi*rand(length(freq),1)-pi; % phases of sum of sines
target2 = sin(freq*2*pi*(0:delt:T-delt) + repmat(phases,1,nstep));
target = sum(target2,1)'; % sum of sines target to track
hand = zeros(nstep,z);
hand(1,:) = -2.5; %initial position of the hand
    
n = 100;
P = zeros(order,order,n);
P(1:order,1:order,1) = rand(order); % use random values for first iteration of P

for i = 2:n
    P(:,:,i) = A'*P(:,:,i-1)*A - (A'*P(:,:,i-1)*B)*inv(R + B'*P(:,:,i-1)*B)*(B'*P(:,:,i-1)*A) + Q;
end
    
L = inv(R + B'*P(:,:,n)*B)*(B'*P(:,:,n)*A); % feedback control law

xt = zeros(order,nstep);
xt(1,1) = -2.5 - target(1); % initialize state variables
hand(1,:) = -2.5; % absolute hand position
u = zeros(nstep,1); % movement commands

for j = 1:z
    for i = 2:nstep
        u(i) = -L*xt(:,i-1);
        xt(:,i) = A*xt(:,i-1) + B*u(i);
        
        hand(i,z) = hand(i-1,z) + (xt(1,i) - xt(1,i-1)); % compute absolute hand position
        xt(1,i) = hand(i,z) - target(i); % adjust xt position according to sum of sines target motion
    end
end

% compute fourier transforms
e = 2/delt; % figure out the number of time steps to throw away

hand = hand((e+1)-lag:(21*e)-lag,:); % time shift hand signal by 'lag'
target = target((e+1):(21*e));
hand_avg = mean(hand,1);
target_avg = mean(target,1);

input_fft = fft(target - repmat(target_avg,[nstep2 1]));
output_fft = fft(hand - repmat(hand_avg,[nstep2 1]));

idx = find(abs(input_fft(:,1))>50); % find the indices of the peaks in the fourier spectrum
idx = idx(1:length(idx)/2);
ratio = output_fft(idx,:)./input_fft(idx,:); % take the complex ratio of output/input
amp = abs(ratio); % magnitude
phase = unwrap(angle(ratio)); % phase

Nblock = size(data.rot.avg.x_x.fft,1);
Nfreq = size(data.rot.avg.x_x.fft,2);

% jitter the frequencies for plotting
a = reshape(datasample([1 -1],Nblock*Nfreq),[Nfreq Nblock]);
x = rand(Nfreq,Nblock).*a*0.01;
y = rand(Nfreq,Nblock).*a*0.01;
scale = repmat(1:Nfreq,[Nblock,1])';
x = x.*scale;
y = y.*scale;
freqs_x_jit = repmat(freqs_x,[Nblock,1])';
freqs_y_jit = repmat(freqs_y,[Nblock,1])';
freqs_x_jit = freqs_x_jit + x;
freqs_y_jit = freqs_y_jit + y;

gains = [0.5 0.75 1];
Yt = 20*log10(gains);
Ytlab = num2cell(gains);

leg = {'Simulation','Empirical'};

%% plot trajectory of hand and target
figure;
plot(t,target,t,hand,'LineWidth',1.5);
%title('Simulation trajectory');
leg2 = {'Target','Hand'};
legend(leg2,'FontSize',15);
xlabel('Time','FontSize',15);
ylabel('Position','FontSize',15);
xlim([0 10]);
% legend('Target','.0001x effort','.01x','Default','100x','10000x');
% legend('Target','10 ms delay','100 ms','200 ms (default)','400 ms','600 ms');

%% plot phasors
figure;
polarplot(ratio,'LineWidth',1.5); hold on;
polarplot(data.rot.avg.x_x.fft(1,:),'-o','LineWidth',1.5);
% polarplot(data.rot.avg.x_x.fft(2,:),'-o','LineWidth',1.5);
% polarplot(data.rot.avg.x_x.fft(5,:),'-o','LineWidth',1.5);
legend(leg,'FontSize',15)

% legend('Inf Horizon','X -> X data','Y -> Y data');
% legend('.0001x effort','.01x','Default','100x','10000x','X -> X data','Y -> Y data');
% legend('10 ms delay','100 ms','200 ms (default)','400 ms','600 ms','X -> X data','Y -> Y data');

%% plot Bode plots
figure; subplot(1,2,1);
semilogx(freq,20*log10(amp),'LineWidth',1.5); hold on;
errorbar(freqs_x_jit(:,1),20*log10(data.rot.avg.x_x.d.amplitude(1,:)),data.rot.avg.x_x.d.amp_err(1,:),'-o','MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',1.5);
% errorbar(freqs_x_jit(:,2),20*log10(data.rot.avg.x_x.d.amplitude(2,:)),data.rot.avg.x_x.d.amp_err(2,:),'-o','MarkerFaceColor',[0.9290 0.6940 0.1250],'LineWidth',3);
% errorbar(freqs_x_jit(:,5),20*log10(data.rot.avg.x_x.d.amplitude(5,:)),data.rot.avg.x_x.d.amp_err(5,:),'-o','MarkerFaceColor',[0.4940 0.1840 0.5560],'LineWidth',3);
% errorbar(freqs_y,20*log10(data.rot.avg.y_y.d.amplitude(1,:)),data.rot.avg.y_y.d.amp_err(1,:),'-o');
set(gca, 'LineWidth', 1, 'box', 'off','Ytick',Yt,'Yticklabel',Ytlab);
title('Gain','FontSize',15);
xlabel('Frequency (Hz)','FontSize',12);
ylabel('Gain (cm/cm)','FontSize',12);
xlim([0 4]);
axis([0.05 2.5 -10 2]);
legend(leg,'FontSize',12,'Location','Southwest','FontSize',12);
pbaspect([1 1 1]);

subplot(1,2,2);
semilogx(freq,unwrap(phase)*180/pi,'LineWidth',2); hold on;
errorbar(freqs_x_jit(:,1),unwrap(data.rot.avg.x_x.d.phase(1,:))*180/pi,data.rot.avg.x_x.d.phase_err(1,:),'-o','MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',2);
% errorbar(freqs_x_jit(:,2),unwrap(data.rot.avg.x_x.d.phase(2,:))*180/pi,data.rot.avg.x_x.d.phase_err(2,:),'-o','MarkerFaceColor',[0.9290 0.6940 0.1250],'LineWidth',3);
% errorbar(freqs_x_jit(:,5),unwrap(data.rot.avg.x_x.d.phase(5,:))*180/pi,data.rot.avg.x_x.d.phase_err(5,:),'-o','MarkerFaceColor',[0.4940 0.1840 0.5560],'LineWidth',3);
% errorbar(freqs_y,unwrap(data.rot.avg.y_y.d.phase(1,:))*180/pi,data.rot.avg.y_y.d.phase_err(1,:),'-o');
set(gca, 'LineWidth',1,'box','off');
title('Phase','FontSize',15);
xlabel('Frequency (Hz)','FontSize',12);
ylabel('Phase (degrees)','FontSize',12);
xlim([0 4]);
axis([0.05 2.5 -350 2]);
pbaspect([1 1 1]);

% legend('Inf horizon','X -> X data','Y -> Y data','Location','southwest');
% legend('.0001x effort','.01x','Default','100x','10000x','X -> X data','Y -> Y data','Location','southwest');
% legend('10 ms delay','100 ms','200 ms (default)','400 ms','600 ms','X -> X data','Y -> Y data','Location','southwest');

%% plot complex tracking error
sim_error = NaN(length(freq),z);
xx_error = NaN(Nblock,Nfreq);
yy_error = NaN(Nblock,Nfreq);

for i = 1:numel(ratio)
    sim_error(i) = norm(1-ratio(i));
end

for i = 1:numel(data.rot.avg.x_x.fft)
    xx_error(i) = norm(1-data.rot.avg.x_x.fft(i));
    yy_error(i) = norm(1-data.rot.avg.y_y.fft(i));
end

figure;
semilogx(freq,sim_error,'LineWidth',3); hold on;
plot(freqs_x,xx_error(1,:),'-o','MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',3);
plot(freqs_x,xx_error(2,:),'-o','MarkerFaceColor',[0.9290 0.6940 0.1250],'LineWidth',3);
plot(freqs_x,xx_error(5,:),'-o','MarkerFaceColor',[0.4940 0.1840 0.5560],'LineWidth',3);
% plot(freqs_y,yy_error,'-o');
title('Tracking Error','FontSize',15);
xlabel('Frequency (Hz)','FontSize',12);
ylabel('Complex Tracking Error','FontSize',12);
axis([0.05 2.5 0 2]);
legend(leg,'FontSize',12,'Location','Northwest');
% legend('.0001x effort','.01x','Default','100x','10000x','X -> X data','Y -> Y data','Location','Northwest');
% legend('10 ms delay','100 ms','200 ms (default)','400 ms','600 ms','X -> X data','Y -> Y data','Location','Northwest');

%% plot response lag
delay = phase./(freq*2*pi);
err = data.rot.avg.x_x.d.phase_err./(360*repmat(freqs_x,[6 1]))*1000;

figure;
plot(freq,-delay*1000,'LineWidth',2); hold on;
errorbar(freqs_x_jit(:,1),-data.rot.avg.x_x.d.delay(1,:)*1000,err(1,:),'-o','MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',2);
% errorbar(freqs_x_jit(:,2),-data.rot.avg.x_x.d.delay(2,:)*1000,err(2,:),'-o','MarkerFaceColor',[0.9290 0.6940 0.1250],'LineWidth',3);
% errorbar(freqs_x_jit(:,5),-data.rot.avg.x_x.d.delay(5,:)*1000,err(5,:),'-o','MarkerFaceColor',[0.4940 0.1840 0.5560],'LineWidth',3);
title('Response Lag','FontSize',15);
xlabel('Frequency (Hz)','FontSize',12);
ylabel('Response Lag (ms)','FontSize',12);
set(gca,'Xscale','log');
axis([0.05 2.5 300 600]);
grid on;
legend(leg,'FontSize',12);
% legend('.0001x effort','.01x','Default','100x','10000x','X -> X data','Y -> Y data');
% legend('10 ms delay','100 ms','200 ms (default)','400 ms','600 ms','X -> X data','Y -> Y data');
