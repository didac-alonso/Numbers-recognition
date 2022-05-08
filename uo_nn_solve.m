%{
Input parameters:

    num_target : digits a identifica.
    tr_freq : freqüencia dels digits en el training set.
    tr_seed : seed per la generació del training set.
    tr_p : tamany del training set.
    te_seed : seed per la generació del test set.
    te_q : tamany del test set.
    la : coeficient de regularització.
    epsG : tolerancia en el òptim, és a dir gL(wo)< epsG.
    kmax : nombre màxim d'iteracions.
    ils : line search usada (1 per exacta, 2 per BLS, 3 per uo_BLSNW32)
    ialmax :  formula pel step lenght màxim (1 or 2).
    kmaxBLS : nombre màxim d'iteracions pel uo_BLSNW32.
    epsal : mínima diferència entre alphas pel uo_BLSNW32
    c1,c2 : paràmetres per les WC
    isd : algorisme d'optimització a usar:
            isd = 1 -> GM
            isd = 3 -> BFGS
            isd = 7 -> SGM
    sg_al0 : \alpha^{SG}_0.
    sg_be : \beta^{SG}.
    sg_ga : \gamma^{SG}.
    sg_emax : e^{SGÇ_{max}.
    sg_ebest : e^{SG}_{best}.
    sg_seed : seed per les permutacions aleatòries del SG.
    icg : if 1 : CGM-FR; if 2, CGM-PR+      (useless in this project).
    irc : re-starting condition for the CGM (useless in this project).
    nu : parameter of the RC2 for the CGM  (useless in this project).

Output parameters:
    Xtr : X^{TR}.
    ytr : y^{TR}.
    wo : w^*.
    fo : {\tilde L}^*.
    tr_acc : Accuracy^{TR}.
    Xte : X^{TE}.
    yte : y^{TE}.
    te_acc : Accuracy^{TE}.
    niter : nombre d'iteracions.
    tex : temps d'execució total (Elapsed time).

%}
function [Xtr,ytr,wo,fo,tr_acc,Xte,yte,te_acc,niter,tex]=uo_nn_solve(num_target,tr_freq,...
    tr_seed,tr_p,te_seed,te_q,la,epsG,kmax,ils,ialmax,kmaxBLS,...
    epsal,c1,c2,isd,sg_al0,sg_be,sg_ga,sg_emax,sg_ebest,sg_seed,icg,irc,nu)    
    
    tic
    % Generació dels datasets:
    [Xtr,ytr] = uo_nn_dataset(tr_seed, tr_p, num_target, tr_freq);
    [Xte,yte] = uo_nn_dataset(te_seed, te_q, num_target, 0);

    
    % Definició les funcions:
    sig = @(X) 1./(1+exp(-X));
    y = @(X,w) sig(w'*sig(X));
    L = @(w) (norm(y(Xtr,w)-ytr)^2)/size(ytr,2)+ (la*norm(w)^2)/2;
    
    % Punt inicial per començar la cerca
    w1 = zeros(35,1);
    
    % Mateixa estructura que el solver
    if isd == 1
        % Inicialització derivada de L
        gL = @(w) (2*sig(Xtr)*((y(Xtr,w)-ytr).*y(Xtr,w).*(1-y(Xtr,w)))')/size(ytr,2)+la*w;
        
        [wk, dk, ak, iWk, niter] = GM(w1, L, gL, 0.001, 1, c1, c2, epsG, kmax, ialmax,kmaxBLS,epsal, ils, 0);
    elseif isd == 3
        % Inicialització derivada de L
        gL = @(w) (2*sig(Xtr)*((y(Xtr,w)-ytr).*y(Xtr,w).*(1-y(Xtr,w)))')/size(ytr,2)+la*w;
        
        [wk, dk, alk, Hk, iWk, niter] = BFGS(w1, L, gL, 0.001, 1, c1, c2, epsG, kmax, ialmax, kmaxBLS, epsal, ils, 0);
    elseif isd == 7
        % Inicialitció L sobre el set de training
        Lte = @(w) (norm(y(Xte,w)-yte)^2)/size(yte,2)+ (la*norm(w)^2)/2;
        % Inicialització de la derivada de L, en aquest cas X i yt són 
        % variables independents
        gL = @(w, X, yt) (2*sig(X)*((y(X,w)-yt).*y(X,w).*(1-y(X,w)))')/size(yt,2)+la*w;

        [wk, niter] = SGM(w1,Xtr, ytr, gL,Lte, sg_al0, sg_be, sg_ga, sg_emax, sg_ebest, sg_seed);
    end
    
    wo = wk(:,end);
    fo = L(wo);
    
    tr_acc = 100*mean(round(y(Xtr, wo))==ytr);

    te_acc = 100*mean(round(y(Xte, wo))==yte);
    tex = toc;
end