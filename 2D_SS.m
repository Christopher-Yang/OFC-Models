clear all;
rng('default');

% adjust time step 
delt = 0.001; % time step length in secs
nstep = ceil(22/delt); % number of time steps
start = [0 0]; % starting position of the hand

% parameters for A and B matrices
t1 = 0.224;
t2 = 0.013;
t3 = 0.004;
k = 0; %1 to include spring
b = t1 + t2;
m = t1*t2;
r = t3;

% sum-of-sines target movement
freq = (0.05:0.05:1)'; % frequencies used in the simulation
phase = 2*pi*rand(length(freq),2)-pi; % phases of sum of sines
target(1,:,:) = sin(freq*2*pi*(0:delt:nstep*delt-delt) + repmat(phase(:,1),1,nstep));
target(2,:,:) = sin(freq*2*pi*(0:delt:nstep*delt-delt) + repmat(phase(:,2),1,nstep));
target = squeeze(sum(target,2));

% generate A, B, and C matrices in discrete time formulation
A = [0 1 0
    -k/m -b/m 1/m
    0 0 -1/r];
B = [0 0 1/r]';
C = [1 0 0];
H = tf(ss(A,B,C,0));
delay = [0 0.5];
ang = [0 0];
name = {'0 s',[num2str(delay(2)),' s']};

for j = 1:length(delay)
    % perform pade approximation
    s = tf('s');
    G = exp(-delay(j)*s);
    G = pade(G,1);
    J = H*G;
    [num,den] = tfdata(J);
    m = length(num{1});

    if delay(j) == 0
        A = zeros(m-1);
        for i = 2:m-1
            A(i-1,i) = 1;
            A(end,i) = -den{1}(end-i+1);
        end
        B = [0 0 num{1}(end)]';
    else
        A = [-den{1}(2) 1 0 0
            -den{1}(3) 0 1 0
            -den{1}(4) 0 0 1
            0 0 0 0];
        B = [0 0 num{1}(4) num{1}(5)]';
        C = [1 0 0 0];
    end

    z = size(A,1)+1;
    Ad = expm(A*delt);
    Ad = blkdiag(Ad, 1);
    Ad = [Ad zeros(z)
         zeros(z) Ad];

    Bd = delt*B;
    Bd = [Bd;0];
    Bd = [Bd zeros(z,1)
         zeros(z,1) Bd];

    C = rotz(ang(j));
    C = C(1:2,1:2);
    C = [diag(repelem(C(1,1),z)) diag(repelem(C(1,2),z))
        diag(repelem(C(2,1),z)) diag(repelem(C(2,2),z))];
    
    % values for Q and R taken from Qian infinite horizon model
    R = [0.1 0
        0 0.1]; % effort cost- default is 0.0001
    Q2 = zeros(z);
    Q2(1,1) = 1;
    Q2(end,end) = 1;
    Q2(1,end) = -1;
    Q2(end,1) = -1;
    Q = [Q2 zeros(z)
        zeros(z) Q2]*0.1;
    
    % state vector
    order = size(Ad,1); % order of the system
    x = zeros(order,nstep);
    x(1,1) = start(1); % x position
    x(end/2+1,1) = start(2); % y position
    x(end/2,1) = target(1,1); % starting target x position
    x(end,1) = target(2,1); % starting target y position
    xhat = x; % predicted state vector
    u = zeros(size(Bd,2),nstep); % movement commands

    % calculate feedback gain
    n = 100000;
    P = zeros(order,order,n);
    P(:,:,1) = rand(order);
    for i = 2:n
        P(:,:,i) = Ad'*P(:,:,i-1)*Ad - (Ad'*P(:,:,i-1)*Bd)*inv(R + Bd'*P(:,:,i-1)*Bd)*(Bd'*P(:,:,i-1)*Ad) + Q;
    end
    L = inv(R + Bd'*P(:,:,i)*Bd)*(Bd'*P(:,:,i)*Ad);

    % simulate trajectory
    for i = 2:nstep
        u(:,i) = -L*inv(C)*x(:,i-1);
        x(:,i) = Ad*x(:,i-1) + Bd*u(:,i);
    %     xhat(:,i) = C*xhat(:,i-1) + Bd*u(:,i-1);

        % set target location
        x(end/2,i) = target(1,i);
        x(end,i) = target(2,i);
        xhat(end/2,i) = target(1,i);
        xhat(end,i) = target(2,i);
    end

    figure(1)
    subplot(1,2,j)
    plot(x(1,:),x(end/2+1,:))
    hold on
    title(['Delay = ',num2str(delay(j)),' s']);
    plot(x(end/2,2:end),x(end,2:end))
    legend({'Hand','Target'})
    pbaspect([1 1 1]);

    % compute fourier transforms
    e = 2/delt; % figure out the number of time steps to throw away

    hand = x([1 end/2+1],(e+1):end)'; % time shift hand signal by 'lag'
    targ = x([end/2 end],(e+1):end)';
    hand_avg = mean(hand,1);
    targ_avg = mean(targ,1);

    fft_in = fft(targ - repmat(targ_avg,[nstep-e 1]));
    fft_out = fft(hand - repmat(hand_avg,[nstep-e 1]));
    [inSort,i] = sort(fft_in(:,1),'descend');
    i = i(1:2*length(freq),:);
    i = sort(i(i<(nstep-e)/2));

    ratio = fft_out(i,:)./fft_in(i,:);
    amp = abs(ratio); % magnitude
    phase = unwrap(angle(ratio)); % phase
    
    figure(2)
    subplot(2,1,1)
    semilogx(freq, amp(:,1))
    hold on
    ylabel('Gain (cm/cm)')
    ylim([0.2 1.5]);
%     title(['Delay = ',num2str(delay(j)),' s']);
    
    subplot(2,1,2)
    semilogx(freq,phase(:,1)*180/pi)
    hold on
    ylabel('Phase (degrees)')
    xlabel('Frequency (Hz)')
%     ylim([-200 0]);
end
figure(2)
legend(name,'Location','Southwest')