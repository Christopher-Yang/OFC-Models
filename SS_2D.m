clear all;
load dat.mat 
load angle.mat
rng(2);

% adjust high level parameters
delt = 0.001; % time step length in secs
start = [0 0]; % starting position of the hand
nReps = 2; % number of base periods to track
Nsubj = 1; % number of subjects to simulate (simulates different phases)
delay = [0.3 0.7]; % amount of delay in seconds
ang = [0 40]; % angle to rotate

% good values
% G = 0.5, I = 0.05, tau = 0.066, delay = 0.25, Q2(1,1) = 1.3 
% R = [0.0015 0; 0 0.0015]

% Single joint reaching movements:
G = 0.14;        % Viscous Constant: Ns/m; originally 0.14
I = 0.1;         % Inertia Kgm2; originally 0.1
tau = 0.066;    % Muscle time constant, s; originally 0.066

% sum-of-sines target movement
freq(:,1) = 0.02:0.02:7; % frequencies used in the simulation
freq(:,2) = 0.03:0.02:7.01;
% freq(:,1) = .1:.1:3;
% freq(:,2) = .15:.1:3.05;
amp = 0.015./freq;
[a1, k1] = min(abs(freq - 0.63));
[a2, k2] = min(a1);
amp(1:k1(k2),:) = amp(k1(k2),k2);
Ts = 1/(min(freq(:))/2);
nstep = ceil(nReps*Ts/delt); % number of time steps

gain = NaN(size(freq,1),4,length(ang),Nsubj);
phase = NaN(size(freq,1),4,length(ang),Nsubj);

for k = 1:Nsubj
    for j = 1:length(ang)
        % create state space model in discrete time
        A = [0 1 0 0
            0 -G/I 4/(delay(j)*I) -1/I
            0 0 -2/delay(j) 1
            0 0 0 -1/tau];
        B = [0 0 0 1/tau]';
        
        z = size(A,1)+1;
        Ad = expm(A*delt);
        Ad = blkdiag(Ad, 1);
        Ad = [Ad zeros(z)
            zeros(z) Ad];
        
        Bd = delt*B;
        Bd = [Bd;0];
        Bd = [Bd zeros(z,1)
            zeros(z,1) Bd];
        
        % accuracy and effort costs
        R = [0.0001 0
            0 0.0001]; % effort cost- default is 0.0001
        if ang(j) == 0
            Q2 = [1.2 0 0 0 -1
                0 0.02 0 0 0
                0 0 0 0 0
                0 0 0 0 0
                -1 0 0 0 1];
            Q = [Q2 zeros(z)
                zeros(z) Q2];
        else
            Q2 = [2.5 0 0 0 -1
                0 0.2 0 0 0
                0 0 0 0 0
                0 0 0 0 0
                -1 0 0 0 1];
            Q = [Q2 zeros(z)
                zeros(z) Q2];
        end
        
        ph = 2*pi*rand(length(freq(:,1)),2)-pi; % phases of sum of sines
        target2(1,:,:) = repmat(amp(:,1),[1 nstep]).*sin(freq(:,1)*2*pi*(0:delt:nstep*delt-delt) + repmat(ph(:,1),[1 nstep]));
        target2(2,:,:) = repmat(amp(:,2),[1 nstep]).*sin(freq(:,2)*2*pi*(0:delt:nstep*delt-delt) + repmat(ph(:,2),[1 nstep]));
        target = squeeze(sum(target2,2));
        
        % rotation matrix
        C = rotz(ang(j));
        C = C(1:2,1:2);
%         C = [diag(repelem(C(1,1),5)) diag(repelem(C(1,2),5))
%             diag(repelem(C(2,1),5)) diag(repelem(C(2,2),5))];
        
        % state vector
        order = size(Ad,1); % order of the system
        x = zeros(order,nstep);
        x(1,1) = start(1); % x position
        x(end/2+1,1) = start(2); % y position
        x(end/2,1) = target(1,1); % starting target x position
        x(end,1) = target(2,1); % starting target y position
        u = zeros(size(Bd,2),nstep); % movement commands
        
        % calculate feedback gain
        n = 5000;
        P = zeros(order,order,n);
        P(:,:,1) = rand(order);
        for i = 2:n
            P(:,:,i) = Ad'*P(:,:,i-1)*Ad - (Ad'*P(:,:,i-1)*Bd)*inv(R + Bd'*P(:,:,i-1)*Bd)*(Bd'*P(:,:,i-1)*Ad) + Q;
        end
        L = inv(R + Bd'*P(:,:,i)*Bd)*(Bd'*P(:,:,i)*Ad)*.2;
        
        % simulate trajectory
        for i = 2:nstep
            if ang(j) == 0
                u(:,i) = -L*x(:,i-1);
            else
                u(:,i) = -C*L*x(:,i-1);
%                 OR
%                 f = rand;
%                 if f < 0.3
%                     u(:,i) = -L*x(:,i-1);
%                 else
%                     u(:,i) = -C*L*x(:,i-1);
%                 end
            end
            x(:,i) = Ad*x(:,i-1) + Bd*u(:,i);
            
            % set target location
            x(end/2,i) = target(1,i);
            x(end,i) = target(2,i);
        end
        
%         figure(1)
%         subplot(1,2,j)
%         plot(x(1,:),x(end/2+1,:))
%         hold on
%         plot(x(end/2,2:end),x(end,2:end))
%         legend({'Hand','Target'})
%         pbaspect([1 1 1])
%         axis([-0.5 0.5 -0.5 0.5])
        
        % compute fourier transforms
        e = round(Ts/delt); % figure out the number of time steps to throw away
        
        hand = x([1 end/2+1],(e+1):end)';
        targ = x([end/2 end],(e+1):end)';
        hand_avg = mean(hand,1);
        targ_avg = mean(targ,1);
        
        fft_in = fft(targ - repmat(targ_avg,[nstep-e 1]));
        fft_out = fft(hand - repmat(hand_avg,[nstep-e 1]));
        [inSort,ix] = sort(fft_in(:,1),'descend'); %ix: x -> x
        [inSort,iy] = sort(fft_in(:,2),'descend'); %iy: y -> y
        ix = ix(1:2*length(freq(:,1)),:);
        iy = iy(1:2*length(freq(:,2)),:);
        ix = sort(ix(ix<(nstep-e)/2));
        iy = sort(iy(iy<(nstep-e)/2));
        
        ratio(:,1) = fft_out(ix,1)./fft_in(ix,1); %x target -> x hand
        ratio(:,2) = fft_out(iy,2)./fft_in(iy,2); %y target -> y hand
        ratio(:,3) = fft_out(iy,1)./fft_in(ix,1); %x target -> y hand
        ratio(:,4) = fft_out(ix,2)./fft_in(iy,2); %y target -> x hand
        gain(:,:,j,k) = abs(ratio); % magnitude
        phase(:,:,j,k) = unwrap(angle(ratio)); % phase
    end
end

gain = mean(gain,4);
phase = mean(phase,4);

figure
subplot(2,1,1)
set(gca,'DefaultLineLineWidth',1.5,'XScale','log','TickDir','out','Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5)
hold on
semilogx(freq(:,1),squeeze(gain(:,1,[1 2])))
% semilogx(freq(:,1),squeeze(gain(:,1,1)))
set(gca,'ColorOrderIndex',1)
semilogx(data.rot.avg.x_x.freqs,data.rot.avg.x_x.amplitude([1 2],:),'--')
% semilogx(data.rot.avg.x_x.d.freqs,data.rot.avg.x_x.d.amplitude(1,:),'--')
ylabel('Gain (cm/cm)')
yticks(0:0.5:1)
title(['Rotation: ',num2str(ang(2)),char(0176)])
legend({['Simulation ' num2str(ang(1)) char(176)],['Simulation ' num2str(ang(2)) char(176)],'Data (Baseline)','Data (Early)'});
% legend({'Simulation','Data (Baseline)'})


subplot(2,1,2)
set(gca,'DefaultLineLineWidth',1.5,'XScale','log','TickDir','out','Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5)
hold on
semilogx(freq(:,1),squeeze(phase(:,1,[1 2]))*180/pi)
% semilogx(freq(:,1),squeeze(phase(:,1,1))*180/pi)
set(gca,'ColorOrderIndex',1)
semilogx(data.rot.avg.x_x.freqs,data.rot.avg.x_x.phase([1 2],:)*180/pi,'--')
% semilogx(data.rot.avg.x_x.d.freqs,data.rot.avg.x_x.d.phase(1,:)*180/pi,'--')
ylabel(['Phase (',char(0176),')'])
xlabel('Frequency (Hz)')
yticks(-360:180:0)
ylim([-400 0])

angle(:,1,1) = atan2(gain(:,3,1),gain(:,1,1));
angle(:,2,1) = atan2(gain(:,4,1),gain(:,2,1));
angle(:,1,2) = atan2(gain(:,1,2),gain(:,3,2));
angle(:,2,2) = atan2(gain(:,2,2),gain(:,4,2));
angle = angle.*180/pi;

figure
set(gca,'DefaultLineLineWidth',2,'XScale','log')
hold on
semilogx(freq(:,1),squeeze(angle(:,1,[1 2])))
set(gca,'ColorOrderIndex',1)
semilogx(data.rot.avg.x_x.freqs,output.rot.X([1 2],:)*180/pi,'--');
title(['Rotation: ',num2str(ang(2)),char(0176)])
ylabel(['Compensation angle (',char(0176),')'])
xlabel('Frequency (Hz)')
yticks(0:45:90)
legend({['Simulation ' num2str(ang(1)) char(176)],['Simulation ' num2str(ang(2)) char(176)],'Data (Baseline)','Data (Early)'},'Location','northwest')
legend boxoff