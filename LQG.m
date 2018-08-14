clear all;
load('dat.mat');
rand('seed',123456);

z = 1; % number of simulations
order = 3; % order of the system
delt = 0.01; % time step in secs

t1 = 0.224; % parameters for A and B
t2 = 0.013;
t3 = 0.004;
k = 1;

b = t1 + t2;
m = t1*t2;
r = t3;

% A = [0 1 0 0 0 0; -k/m -b/m 1/m 0 0 0; 0 0 -1/r 0 0 0; 0 0 0 0 1 0; 0 0 0 -k/m -b/m 1/m; 0 0 0 0 0 -1/r];
% A = [0 1 0 0; -k/m -b/m 1/m 0; 0 0 -1/r 0; 0 0 0 1];
A = [0 1 0; -k/m -b/m 1/m; 0 0 -1/r];
% a = [0 1 0; -k/m -b/m 1/m; 0 0 -1/r];
% a = eye(3) + delt*a;
% A = [a zeros(3); eye(3) zeros(3)];
B = [0 0 1/r]';
B = delt*B;
% C = [zeros(3) eye(3)];

T = 42; % total simulation time
T2 = 40; % amount of analysis data
t = 0:delt:T2-delt; % x axis for graphs
nstep = round(T/delt); % number of simulation time steps
nstep2 = round(T2/delt); % number of analysis time steps

freq = (0.05:0.05:4)'; % frequencies used in the simulation
freqs_x = data.rot.avg.x_x.d.freqs; % frequencies of experimental x data
freqs_y = data.rot.avg.y_y.d.freqs; % frequencies of experimental y data
phases = 2*pi*rand(length(freq),1)-pi; % phases of sum of sines
target2 = sin(freq*2*pi*(0:delt:T-delt) + repmat(phases,1,nstep));
target = sum(target2,1)'; % sum of sines target to track
hand = zeros(nstep,z); 
hand(1) = -2.5; %initial position of the hand

% desired: (x-g)^2 + x^2 + u^2

% Q = diag([1 0.01 0]); % accuracy cost- originally diag([1 0.01 0])
% Q = [2 0 0 -1; 0 0 0 0; 0 0 0 0; -1 0 0 1];
Q = diag([1 0 0]);
R = 0.0001; % effort cost- originally 0.0001
% lag = [1 10 20 40 60];
lag = 1;

for j = 1:z
%     if j ~= 1
%         Q = Q*100;
%     end

    n = 100;
    P = zeros(order,order,n);
    P(1:order,1:order,1) = rand(order); % for one dimensional case
%     P(1:order/2,1:order/2,1) = rand(order/2); % for two dimensional case
%     P(order/2+1:end,order/2+1:end,1) = P(1:order/2,1:order/2,1);
    
    for i = 2:n
        P(:,:,i) = A'*P(:,:,i-1)*A - (A'*P(:,:,i-1)*B)*inv(R + B'*P(:,:,i-1)*B)*(B'*P(:,:,i-1)*A) + Q;
%         L(:,i) = inv(R + B'*P(:,:,i)*B)*(B'*P(:,:,i)*A);
    end
    
    L = inv(R + B'*P(:,:,i)*B)*(B'*P(:,:,i)*A);
%     L = [L(4:6) L(1:3)];
%     L = L(:,i)';
    
    xt = zeros(order,nstep);
    xhatt = zeros(order,nstep);
%     xt(1,1:lag) = -2.5;
%     hand(1:lag,j) = -2.5;
    xt(1,1) = -2.5 - target(1);
%     xt(4,1) = -2.5 - target(1);
    xhatt(1,1) = -2.5 - target(1);
%     xhatt(4,1) = -2.5 - target(1);
    u = zeros(nstep,1);
    for i = 2:nstep
        u(i) = -L*xt(:,i-1);
        xt(:,i) = A*xt(:,i-1) + B*u(i);

        hand(i) = hand(i-1) + (xt(1,i) - xt(1,i-1));
        xt(1,i) = hand(i) - target(i);
        
%         dy = C*xt(:,i-1);
%         xhatt(:,i) = 

    end
end

% compute fourier transforms
hand = hand(201:4200,:);
target = target(201:4200);
hand_avg = mean(hand,1);
target_avg = mean(target,1);

input_fft = fft(target - repmat(target_avg,[nstep2 1]));
output_fft = fft(hand - repmat(hand_avg,[nstep2 1]));

idx = find(abs(input_fft(:,1))>50);
idx = idx(1:length(idx)/2);
ratio = output_fft(idx,:)./input_fft(idx,:);
amp = abs(ratio);
phase = unwrap(angle(ratio));

%% plot trajectory of hand and target
figure;
plot(t,target,t,hand);
legend('Target','Hand');
xlabel('Time'); 
ylabel('Position');
% xlim([0 10]);
% legend('Target','.0001x effort','.01x','Default','100x','10000x');
% legend('Target','10 ms delay','100 ms','200 ms (default)','400 ms','600 ms');

%% plot phasors
figure;
polarplot(ratio); hold on;
polarplot(data.rot.avg.x_x.fft(1,:),'-o');
polarplot(data.rot.avg.y_y.fft(1,:),'-o');

% legend('Inf Horizon','X -> X data','Y -> Y data');
% legend('.0001x effort','.01x','Default','100x','10000x','X -> X data','Y -> Y data');
% legend('10 ms delay','100 ms','200 ms (default)','400 ms','600 ms','X -> X data','Y -> Y data');

%% plot Bode plots
figure; subplot(2,1,1);
plot(freq,20*log10(amp)); hold on;
errorbar(freqs_x,20*log10(data.rot.avg.x_x.d.amplitude(1,:)),data.rot.avg.x_x.d.amp_err(1,:),'-o');
errorbar(freqs_y,20*log10(data.rot.avg.y_y.d.amplitude(1,:)),data.rot.avg.y_y.d.amp_err(1,:),'-o');
set(gca,'Xscale','log');
title('Bode Plots');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([0 4]);

subplot(2,1,2);
plot(freq,unwrap(phase)*180/pi); hold on;
errorbar(freqs_x,unwrap(data.rot.avg.x_x.d.phase(1,:))*180/pi,data.rot.avg.x_x.d.phase_err(1,:),'-o');
errorbar(freqs_y,unwrap(data.rot.avg.y_y.d.phase(1,:))*180/pi,data.rot.avg.y_y.d.phase_err(1,:),'-o');
set(gca,'Xscale','log');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
xlim([0 4]);
% legend('Inf horizon','X -> X data','Y -> Y data','Location','southwest');
% legend('.0001x effort','.01x','Default','100x','10000x','X -> X data','Y -> Y data','Location','southwest');
% legend('10 ms delay','100 ms','200 ms (default)','400 ms','600 ms','X -> X data','Y -> Y data','Location','southwest');

%% plot complex tracking error
sim_error = NaN(length(freq),z);
xx_error = NaN(length(freqs_x),1);
yy_error = NaN(length(freqs_y),1);

for i = 1:numel(ratio)
    sim_error(i) = norm(1-ratio(i));
end

for i = 1:length(freqs_x)
    xx_error(i) = norm(1-data.rot.avg.x_x.fft(1,i));
    yy_error(i) = norm(1-data.rot.avg.y_y.fft(1,i));
end

figure;
% semilogx(freq,sim_error,'-o');
semilogx(freq,sim_error,freqs_x,xx_error,'-o',freqs_y,yy_error,'-o');
title('Tracking Error');
xlabel('Frequency (Hz)');
ylabel('Complex Tracking Error');
% legend('.0001x effort','.01x','Default','100x','10000x','X -> X data','Y -> Y data','Location','Northwest');
% legend('10 ms delay','100 ms','200 ms (default)','400 ms','600 ms','X -> X data','Y -> Y data','Location','Northwest');
% legend('Inf Horizon','X -> X data','Y -> Y data');

%% response time
delay = phase./(freq*2*pi);

figure;
plot(freq,-delay*1000,freqs_x,-data.rot.avg.x_x.d.delay(1,:)*1000,'-o',freqs_y,-data.rot.avg.y_y.d.delay(1,:)*1000,'-o');
title('Response Time');
xlabel('Frequency (Hz)');
ylabel('Response Time (ms)');
% legend('.0001x effort','.01x','Default','100x','10000x','X -> X data','Y -> Y data');
% legend('10 ms delay','100 ms','200 ms (default)','400 ms','600 ms','X -> X data','Y -> Y data');
