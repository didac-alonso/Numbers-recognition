function [Xtr,ytr,wo,fo,tr_acc,Xte,yte,te_acc,niter,tex]=uo_nn_solve(num_target,tr_freq,...
    tr_seed,tr_p,te_seed,te_q,la,epsG,kmax,ils,ialmax,kmaxBLS,...
    epsal,c1,c2,isd,sg_al0,sg_be,sg_ga,sg_emax,sg_ebest,sg_seed,icg,irc,nu)    
     
    tic
    % Generem el dataset de training:
    [Xtr,ytr] = uo_nn_dataset(tr_seed, tr_p, num_target, tr_freq);
    % Definim les funcions:
    sig = @(X) 1./(1+exp(-X));
    y = @(X,w) sig(w'*sig(X));
    L = @(w) (norm(y(Xtr,w)-ytr)^2)/size(ytr,2)+ (la*norm(w)^2)/2;
    [Xte,yte] = uo_nn_dataset(te_seed, te_q, num_target, 0);
    % Punt inicial per començar la cerca
    w1 = zeros(35,1);
    % Mateixa estructura que el solver
    if isd == 1
        gL = @(w) (2*sig(Xtr)*((y(Xtr,w)-ytr).*y(Xtr,w).*(1-y(Xtr,w)))')/size(ytr,2)+la*w;
        
        [wk, dk, ak, iWk, niter] = GM(w1, L, gL, 1, c1, c2, epsG, kmax, ialmax,kmaxBLS,epsal);
    elseif isd == 2
        [wk, dk, alk, betak, iWk, niter] = CGM(w1, L, gL, sig, y, almin, almax, rho, c1, c2, iW, epsG, kmax, icg, irc, nu);
    elseif isd == 3
        gL = @(w) (2*sig(Xtr)*((y(Xtr,w)-ytr).*y(Xtr,w).*(1-y(Xtr,w)))')/size(ytr,2)+la*w;
        
        [wk, dk, alk, Hk, iWk, niter] = BFGS(w1, L, gL, 1, c1, c2, epsG, kmax, ialmax, kmaxBLS, epsal);
    elseif isd == 7
        Lte = @(w) (norm(y(Xte,w)-yte)^2)/size(yte,2)+ (la*norm(w)^2)/2;
        gL = @(w, X, yt) (2*sig(X)*((y(X,w)-yt).*y(X,w).*(1-y(X,w)))')/size(yt,2)+la*w;

        [wk, niter] = SGM(w1,Xtr, ytr, L, gL,Lte, sg_al0, sg_be, sg_ga, sg_emax, sg_ebest, sg_seed);
    end
    % Això faré que el GM no torni un vector directament.
    wo = wk(:,end);
    fo = L(wo);
    tr_acc = 100*mean(round(y(Xtr, wo))==ytr);
    % Generem set de test
    te_acc = 100*mean(round(y(Xte, wo))==yte);
    tex = toc;
end