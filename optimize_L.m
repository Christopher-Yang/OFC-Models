% optimize gain
Linit = [50 10 5];
tau = 1000*[0.324 0.33 0.2];

f_targ = @(L) sim_error(L,tau);

Lopt = fmincon(f_targ,Linit,[],[]);