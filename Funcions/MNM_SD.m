%{
Input:
    - delta: valor mínim que tindrem a la diagonal de la matriu de vaps de
    la hessiana
Output:
    - Hk: hessianes modificades
%}
function [xk, dk, alk, Hk, iWk, it]= MNM_SD(x, f, df, d2f, amin, amax, p, c1, c2, iW, tol,itmax, delta, Q)
    it=1;
    xk = [x];
    dfk = df(x);
    alk = []; dk = []; iWk = []; Hk = [];
    while norm(dfk)>tol & it<=itmax
        [B,Binv, nmod] = SD_mod(d2f(x),delta);
        d = -Binv*dfk;
        [a, iWout] = BLS(x, f, df, d, amin, amax, p, c1, c2, iW, Q);

        x=x+a*d;
        dfk = df(x);

        xk = [xk x];        
        Hk = cat(3, Hk, B);
        dk = [dk d];
        alk = [alk a];
        iWk = [iWk iWout*(~nmod)+4*nmod*(a==1)]; % 4 si usem Newton
        it=it+1;
    end    
end