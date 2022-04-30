clear;
%
% Parameters for dataset generation
%
num_target =[10];
tr_freq    = .5;        
tr_p       = 250;       
te_q       = 250;       
tr_seed    = 123456;    
te_seed    = 789101;    
%
% Parameters for optimization
%
la = 0.0;                                                     % L2 regularization.
epsG = 10^-6; kmax = 10000;                                   % Stopping criterium.
ils=3; ialmax = 2; kmaxBLS=30; epsal=10^-3;c1=0.01; c2=0.45;  % Linesearch.
isd = 1; icg = 2; irc = 2 ; nu = 1.0;                         % Search direction.
sg_seed = 565544; sg_al0 = 2; sg_be = 0.3; sg_ga = 0.01;      % SGM iteration.
sg_emax = kmax; sg_ebest = floor(0.01*sg_emax);               % SGM stopping condition.
%
% Optimization
%
t1=clock;
[Xtr,ytr,wo,fo,tr_acc,Xte,yte,te_acc,niter,tex]=uo_nn_solve(num_target,tr_freq,tr_seed,tr_p,te_seed,te_q,la,epsG,kmax,ils,ialmax,kmaxBLS,epsal,c1,c2,isd,sg_al0,sg_be,sg_ga,sg_emax,sg_ebest,sg_seed,icg,irc,nu);
t2=clock;
fprintf(' wall time = %6.1d s.\n', etime(t2,t1));
%
[Xtr,ytr] = uo_nn_dataset(tr_seed, tr_p, num_target, tr_freq);
% Definim les funcions:
sig = @(X) 1./(1+exp(-X));
y = @(X,w) sig(w'*sig(X));
L = @(w) (norm(y(Xtr,w)-ytr)^2)/size(ytr,2)+ (la*norm(w)^2)/2;
gL = @(w) (2*sig(Xtr)*((y(Xtr,w)-ytr).*y(Xtr,w).*(1-y(Xtr,w)))')/size(ytr,2)+la*w;

% Train
figure(1)
uo_nn_Xyplot(Xtr,ytr,wo)

% Test
[Xte,yte] = uo_nn_dataset(te_seed, te_q, num_target, 0);
figure(2)
uo_nn_Xyplot(Xte,yte,wo)
%% Set 2
% Training data set generation.
clear
la = 0;
epsG= 1.0e-06; kmax= 10000;
ils= 3; ialmax= 2; kmaxBLS= 30; epsal= 1.0e-03;
c1= 0.01; c2= 0.45; isd= 1;
icg = 0; nu = 0; irc = 0;
sg_seed = 565544; sg_al0 = 2; sg_be = 0.3; sg_ga = 0.01;      % SGM iteration.
sg_emax = kmax; sg_ebest = floor(0.01*sg_emax);               % SGM stopping condition.



num_target = [3];
tr_freq = 0.50;
tr_p = 250;
tr_seed = 123456;
% Test data set generation.
te_freq = 0.00;
te_q = 250;
te_seed = 789101;

t1=clock;
[Xtr,ytr,wo,fo,tr_acc,Xte,yte,te_acc,niter,tex]=uo_nn_solve(num_target,tr_freq,tr_seed,tr_p,te_seed,te_q,la,epsG,kmax,ils,ialmax,kmaxBLS,epsal,c1,c2,isd,sg_al0,sg_be,sg_ga,sg_emax,sg_ebest,sg_seed,icg,irc,nu);
t2=clock;
fprintf(' wall time = %6.1d s.\n', etime(t2,t1));

% Train
figure(1)
uo_nn_Xyplot(Xtr,ytr,wo)

% Test
figure(2)
uo_nn_Xyplot(Xte,yte,wo)

%% Exemple 3:
% Training data set generation.
clear
la = 0;
epsG= 1.0e-06; kmax= 10000;
ils= 3; ialmax= 2; kmaxBLS= 30; epsal= 1.0e-03;
c1= 0.01; c2= 0.45; isd= 1;
num_target = [1 3 5 7];
tr_freq = 0.50; icg = 0;
tr_p = 250; irc = 0; nu = 0;
tr_seed = 123456;
% Test data set generation.
te_freq = 0.00;
te_q = 250;
te_seed = 789101;

sg_seed = 565544; sg_al0 = 2; sg_be = 0.3; sg_ga = 0.01;      % SGM iteration.
sg_emax = kmax; sg_ebest = floor(0.01*sg_emax);               % SGM stopping condition.



t1=clock;
[Xtr,ytr,wo,fo,tr_acc,Xte,yte,te_acc,niter,tex]=uo_nn_solve(num_target,tr_freq,tr_seed,tr_p,te_seed,te_q,la,epsG,kmax,ils,ialmax,kmaxBLS,epsal,c1,c2,isd,sg_al0,sg_be,sg_ga,sg_emax,sg_ebest,sg_seed,icg,irc,nu);
t2=clock;
fprintf(' wall time = %6.1d s.\n', etime(t2,t1));
% Train
figure(1)
uo_nn_Xyplot(Xtr,ytr,wo)

% Test
figure(2)
uo_nn_Xyplot(Xte,yte,wo)

%%
clear;
%
% Parameters for dataset generation
%
num_target =[1];
tr_freq    = .5;        
tr_p       = 250;       
te_q       = 250;       
tr_seed    = 123456;    
te_seed    = 789101;    
%
% Parameters for optimization
%
la = 0.0;                                                     % L2 regularization.
epsG = 10^-6; kmax = 10000;                                   % Stopping criterium.
ils=3; ialmax = 2; kmaxBLS=30; epsal=10^-3;c1=0.01; c2=0.45;  % Linesearch.
isd = 3; icg = 2; irc = 2 ; nu = 1.0;                         % Search direction.
sg_seed = 565544; sg_al0 = 2; sg_be = 0.3; sg_ga = 0.01;      % SGM iteration.
sg_emax = kmax; sg_ebest = floor(0.01*sg_emax);               % SGM stopping condition.
%
% Optimization
%
t1=clock;
[Xtr,ytr,wo,fo,tr_acc,Xte,yte,te_acc,niter,tex]=uo_nn_solve(num_target,tr_freq,tr_seed,tr_p,te_seed,te_q,la,epsG,kmax,ils,ialmax,kmaxBLS,epsal,c1,c2,isd,sg_al0,sg_be,sg_ga,sg_emax,sg_ebest,sg_seed,icg,irc,nu);
t2=clock;
fprintf(' wall time = %6.1d s.\n', etime(t2,t1));
%
[Xtr,ytr] = uo_nn_dataset(tr_seed, tr_p, num_target, tr_freq);
% Definim les funcions:
sig = @(X) 1./(1+exp(-X));
y = @(X,w) sig(w'*sig(X));
L = @(w) (norm(y(Xtr,w)-ytr)^2)/size(ytr,2)+ (la*norm(w)^2)/2;
gL = @(w) (2*sig(Xtr)*((y(Xtr,w)-ytr).*y(Xtr,w).*(1-y(Xtr,w)))')/size(ytr,2)+la*w;

% Train
figure(1)
uo_nn_Xyplot(Xtr,ytr,wo)

% Test
[Xte,yte] = uo_nn_dataset(te_seed, te_q, num_target, 0);
figure(2)
uo_nn_Xyplot(Xte,yte,wo)
%% Set 2
% Training data set generation.
clear
la = 0;
epsG= 1.0e-06; kmax= 10000;
ils= 3; ialmax= 1; kmaxBLS= 30; epsal= 1.0e-03;
c1= 0.01; c2= 0.45; isd= 7;
icg = 0; nu = 0; irc = 0;
sg_seed = 565544; sg_al0 = 2; sg_be = 0.3; sg_ga = 0.01;      % SGM iteration.
sg_emax = kmax; sg_ebest = floor(0.01*sg_emax);               % SGM stopping condition.



num_target = [3];
tr_freq = 0.50;
tr_p = 250;
tr_seed = 123456;
% Test data set generation.
te_freq = 0.00;
te_q = 250;
te_seed = 789101;

t1=clock;
[Xtr,ytr,wo,fo,tr_acc,Xte,yte,te_acc,niter,tex]=uo_nn_solve(num_target,tr_freq,tr_seed,tr_p,te_seed,te_q,la,epsG,kmax,ils,ialmax,kmaxBLS,epsal,c1,c2,isd,sg_al0,sg_be,sg_ga,sg_emax,sg_ebest,sg_seed,icg,irc,nu);
t2=clock;
fprintf(' wall time = %6.1d s.\n', etime(t2,t1));

% Train
figure(1)
uo_nn_Xyplot(Xtr,ytr,wo)

% Test
figure(2)
uo_nn_Xyplot(Xte,yte,wo)

%% Exemple 3:
% Training data set generation.
clear
la = 0;
epsG= 1.0e-06; kmax= 10000;
ils= 3; ialmax= 1; kmaxBLS= 30; epsal= 1.0e-03;
c1= 0.01; c2= 0.45; isd= 7;
num_target = [1 3 5 7 9];
tr_freq = 0.50; icg = 0;
tr_p = 250; irc = 0; nu = 0;
tr_seed = 123456;
% Test data set generation.
te_freq = 0.00;
te_q = 250;
te_seed = 789101;

sg_seed = 565544; sg_al0 = 2; sg_be = 0.3; sg_ga = 0.01;      % SGM iteration.
sg_emax = kmax; sg_ebest = floor(0.01*sg_emax);               % SGM stopping condition.



t1=clock;
[Xtr,ytr,wo,fo,tr_acc,Xte,yte,te_acc,niter,tex]=uo_nn_solve(num_target,tr_freq,tr_seed,tr_p,te_seed,te_q,la,epsG,kmax,ils,ialmax,kmaxBLS,epsal,c1,c2,isd,sg_al0,sg_be,sg_ga,sg_emax,sg_ebest,sg_seed,icg,irc,nu);
t2=clock;
fprintf(' wall time = %6.1d s.\n', etime(t2,t1));
tr_acc
te_acc
% Train
figure(1)
uo_nn_Xyplot(Xtr,ytr,wo)

% Test
figure(2)
uo_nn_Xyplot(Xte,yte,wo)

%% EXEMPLE SGM:
clear;
%
% Parameters for dataset generation
%
num_target =[3];
tr_freq = .5;
tr_p = 250;
te_q = 250;
tr_seed = 123456;
te_seed = 789101;
%
% Parameters for optimization
%
la = 0.0; % L2 regularization.
epsG = 10^-6; kmax = 10000; % Stopping criterium.
ils=3; ialmax = 2; kmaxBLS=30; epsal=10^-3;c1=0.01; c2=0.45; % Linesearch.
isd = 7; icg = 2; irc = 2 ; nu = 1.0; % Search direction.
sg_seed = 565544; sg_al0 = 2; sg_be = 0.3; sg_ga = 0.01; % SGM iteration.
sg_emax = kmax; sg_ebest = floor(0.01*sg_emax); % SGM stopping condition.
%
% Optimization
%
t1=clock;
[Xtr,ytr,wo,fo,tr_acc,Xte,yte,te_acc,niter,tex]=uo_nn_solve(num_target,tr_freq,tr_seed,tr_p,te_seed,te_q,la, ...
epsG,kmax,ils,ialmax,kmaxBLS,epsal,c1,c2,isd,sg_al0,sg_be,sg_ga,sg_emax,sg_ebest,sg_seed,icg,irc,nu);
t2=clock;
fprintf(' wall time = %6.1d s.\n', etime(t2,t1));
tr_acc
te_acc
niter
% Train
figure(1)
uo_nn_Xyplot(Xtr,ytr,wo)

% Test
figure(2)
uo_nn_Xyplot(Xte,yte,wo)
