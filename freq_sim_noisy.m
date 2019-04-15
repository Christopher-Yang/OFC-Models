function sim = freq_sim_noisy(L,i_seed,tau)

%load('av_data.mat');
rng(1); % set seed for phase randomization

z = 1; % number of simulations
order = 3; % order of the system
delt = 0.005; % time step in secs

% adjust amount of time lag 
lag = 0.27/delt; % compute number of time steps of lag (0.27 = 270 milliseconds)

% values for Q and R taken from Qian infinite horizon model
Q = diag([1 0.01 0]); % accuracy cost- default is [1 0.01 0]
R = 0.0001; % effort cost- default is 0.0001

% parameters for A and B matrices
if(nargin<=2)
    t1 = 0.224;
    t2 = 0.13;
    t3 = 0.04;
else
    t1 = tau(1); t2 = tau(2); t3 = tau(3);
end

k = 0;
b = t1 + t2;
m = t1*t2;
r = t3;

sigma = 10000*sqrt(delt);

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

%freqs_x = data.rot.avg.x_x.d.freqs; % frequencies of experimental x data
%freqs_y = data.rot.avg.y_y.d.freqs; % frequencies of experimental y data

freq = (0.05:0.05:2.5)'; % frequencies used in the simulation
phases = 2*pi*rand(length(freq),1)-pi; % phases of sum of sines
target2 = sin(freq*2*pi*(0:delt:T-delt) + repmat(phases,1,nstep));
target = sum(target2,1)'; % sum of sines target to track
hand = zeros(nstep,z);
hand(1,:) = -2.5; %initial position of the hand
   
%{
n = 100;
P = zeros(order,order,n);
P(1:order,1:order,1) = rand(order); % use random values for first iteration of P

for i = 2:n
    P(:,:,i) = A'*P(:,:,i-1)*A - (A'*P(:,:,i-1)*B)*inv(R + B'*P(:,:,i-1)*B)*(B'*P(:,:,i-1)*A) + Q;
end
%}  
%L = inv(R + B'*P(:,:,n)*B)*(B'*P(:,:,n)*A); % feedback control law

%L(1) = 40;
%L(2) = 4;

%L = [75 7 4];

xt = zeros(order,nstep);
xt(1,1) = -2.5 - target(1); % initialize state variables
hand(1,:) = -2.5; % absolute hand position
u = zeros(nstep,1); % movement commands

rng(i_seed); % set seed for noise during execution

for j = 1:z
    for i = 2:nstep
        u(i) = -L*xt(:,i-1);
        xt(:,i) = A*xt(:,i-1) + B*(u(i)+sigma*randn(1));
        
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

%Nblock = size(data.rot.avg.x_x.fft,1);
%Nfreq = size(data.rot.avg.x_x.fft,2);

% jitter the frequencies for plotting
%a = reshape(datasample([1 -1],Nblock*Nfreq),[Nfreq Nblock]);
%x = rand(Nfreq,Nblock).*a*0.01;
%y = rand(Nfreq,Nblock).*a*0.01;
%scale = repmat(1:Nfreq,[Nblock,1])';
%x = x.*scale;
%y = y.*scale;
%freqs_x_jit = repmat(freqs_x,[Nblock,1])';
%freqs_y_jit = repmat(freqs_y,[Nblock,1])';
%freqs_x_jit = freqs_x_jit + x;
%freqs_y_jit = freqs_y_jit + y;

%gains = [0.5 0.75 1];
%Yt = 20*log10(gains);
%Ytlab = num2cell(gains);

sim.t = t;
sim.output_fft = output_fft;
sim.input_fft = input_fft;
sim.hand = hand;
sim.target = target;
sim.idx = idx;
sim.ratio = ratio;
%sim.freqs_x = freqs_x;
%sim.freqs_y = freqs_y;
sim.freq = freq;
sim.L = L;
