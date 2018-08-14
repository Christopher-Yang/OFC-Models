clear all;
load('dat.mat');
rand('seed',1);

delt = 0.0005; % time step in secs
lag = 0/delt; % compute number of time steps of lag (0.27 = 270 milliseconds)

% values for Q and R taken from Qian infinite horizon model
Q = diag([1 0.1 0.01 0.0001]); % accuracy cost- default is [1 0.01 0]
R = 0.0001; % effort cost- default is 0.0001

t1 = 0.224; 
t2 = 0.013;
t3 = 0.004;
k = 0;
b = t1 + t2;
m = t1*t2;
r = t3;

% generate A and B matrices in discrete time formulation
Aa = [0 1 0; -k/m -b/m 1/m; 0 0 -1/r];
% A = [Aa 0*Aa; 0*Aa Aa];
% B = [0 0 1/r 0 0 0; 0 0 0 0 0 1/r]';
B = [0 0 1/r]';
order = size(B,1); % order of the system
C = [1 0 0];
D = 0;

% calculate transfer function
delay = 0.27;
if delay == 0
    f1 = 1;
    f2 = 0;
else
    f1 = 1;
    f2 = delay/2;
end
syms s
h1 = det([s*eye(3)-Aa -B; C D]); % numerator symbolically
h2 = det(s*eye(3) - Aa); % denominator symbolically
p1 = f1 - f2*s;
p2 = f1 + f2*s;
num = coeffs(p1*h1);
den = coeffs(p2*h2);
num = double(fliplr(num));
den = double(fliplr(den));
num = num*(1/den(1));
den = den*(1/den(1));
den(5) = 0; %% use when k = 0 (no spring constant)
H = tf(num, den);
G = ss(H);
[A, B, C, D] = ssdata(G);
Ad = eye(order) + delt*A;
Bd = delt*B;

[L,P,e] = lqr(G,Q,R); 

T = 42; % total simulation time
T2 = 40; % amount of analysis data
t = 0:delt:T2-delt; % x axis for graphs
nstep = round(T/delt); % number of simulation time steps
nstep2 = round(T2/delt); % number of analysis time steps

target = NaN(nstep,size(R,1));
freq = (0.05:0.05:2.5)'; % frequencies used in the simulation
phases_x = 2*pi*rand(length(freq),1)-pi; % phases of sum of sines
target_x = sin(freq*2*pi*(0:delt:T-delt) + repmat(phases_x,1,nstep));
target(:,1) = sum(target_x,1)'; % sum of sines target to track

% phases_y = 2*pi*rand(length(freq),1)-pi;
% target_y = sin(freq*2*pi*(0:delt:T-delt) + repmat(phases_y,1,nstep));
% target(:,2) = sum(target_y,1)';
%%
hand = zeros(nstep,size(R,1)); 
hand(1,:,:) = -2.5; %initial position of the hand

X = zeros(order,nstep);
Xhat = zeros(order,nstep);
X(1,1) = -2.5 - target(1,1); % initialize state variables; x position
% xt(4,1) = -2.5 - target(1,2); % y position
Xhat(1,1) = -2.5 - target(1,1);
u = zeros(size(R,1),nstep); % movement commands

for i = 2:nstep
    u(:,i) = -L*Xhat(:,i-1);
    X(:,i) = Ad*X(:,i-1) + Bd*u(:,i);
    
    hand(i,1) = hand(i-1,1) + (X(1,i) - X(1,i-1)); % compute absolute hand position
%     hand(i,2,z) = hand(i-1,2,z) + (X(4,i) - X(4,i-1));
    
    X(1,i) = hand(i,1,z) - target(i,1); % adjust xt position according to sum of sines target motion
    X(4,i) = hand(i,2,z) - target(i,2);
end

% compute fourier transforms
% e = 2/delt; % figure out the number of time steps to throw away
% 
% hand = hand((e+1):(21*e),:,:); % time shift hand signal by 'lag'
% target = target((e+1):(21*e),:);
% hand_avg = mean(hand,1);
% target_avg = mean(target,1);
% 
% input_fft = fft(target - repmat(target_avg,[nstep2 1]));
% output_fft = fft(hand - repmat(hand_avg,[nstep2 1]));
% 
% idx = find(abs(input_fft(:,1))>50); % find the indices of the peaks in the fourier spectrum
% idx = idx(1:length(idx)/2);
% ratio = output_fft(idx,:,:)./repmat(input_fft(idx,:),[1,1,length(lag)]); % take the complex ratio of output/input
% amp = abs(ratio); % magnitude
% phase = unwrap(angle(ratio)); % phase
% 
Nblock = size(data.rot.avg.x_x.fft,1);
Nfreq = size(data.rot.avg.x_x.fft,2);

% jitter the frequencies for plotting
% a = reshape(datasample([1 -1],Nblock*Nfreq),[Nfreq Nblock]);
% x = rand(Nfreq,Nblock).*a*0.01;
% y = rand(Nfreq,Nblock).*a*0.01;
% scale = repmat(1:Nfreq,[Nblock,1])';
% x = x.*scale;
% y = y.*scale;
% freqs_x_jit = repmat(freqs_x,[Nblock,1])';
% % freqs_y_jit = repmat(freqs_y,[Nblock,1])';
% freqs_x_jit = freqs_x_jit + x;
% % freqs_y_jit = freqs_y_jit + y;

% leg = {['Sim ',num2str(ang(1)),'^{\circ}'],['Sim ',num2str(ang(2)),'^{\circ}'],['Sim ',num2str(ang(3)),'^{\circ}'],'Baseline','No training','Max training'};

freqs_x = data.rot.avg.x_x.d.freqs;
% o = bodeplot(H); hold on;
% e = getoptions(o);
% e.PhaseMatching = 'on';
% e.PhaseMatchingFreq = 0.05;
% e.PhaseMatchingValue = 0;
% setoptions(o,e);
% plot(freqs_x,180/pi*(data.rot.avg.x_x.d.phase(1,:)),'-ro');
% axis([0.09 2.1 -360 0]);

figure;
step(H)
figure;
pzmap(H)
lsim(H,target(:,1),0:delt:T-delt);
lsim(G,target(:,1),0:delt:T-delt);
