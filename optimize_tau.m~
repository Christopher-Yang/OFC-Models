% optimize gain
L = [50 10 5];
L = [46 35 4.8];
tau_init = 1000*[0.324 0.2 0.1];

f_targ = @(tau) sim_error(L,tau);

Lopt = fmincon(f_targ,tau_init,[],[],[],[],zeros(1,3),1000*ones(1,3));