%{
Usa el mètode CMI (cholesky) per modificar el newton. I convergir a minims.

Output:
    - Hk: conté les Hessianes modificades
    - tauk : vector amb cada tau usada per modificar la hessiana.
        si tau = 0, no la hem modificat
%}
function [xk, dk, alk, Hk, iWk, tauk, it]= MNM_CMI(x, f, df, d2f, amin, amax, p, c1, c2, iW, tol,itmax, Q)
    it=1;
    xk = [x];
    dfk = df(x);
    Hk = []; dk = []; alk = []; iWk = []; tauk = [];
    while norm(dfk)>tol & it<=itmax
        [B,tau] = CMI_mod(d2f(x));
        d = -B\dfk; % B\ = B^-1*
        [a, iWout] = BLS(x, f, df, d, amin, amax, p, c1, c2, iW, Q);
        
        x=x+a*d;
        dfk = df(x);

        xk = [xk x];
        dk = [dk d];
        alk = [alk a];
        iWk = [iWk (tau~=0)*iWout+(tau==0)*4*(a==1)]; % 4 si és pas de Newton
        tauk = [tauk tau];
        Hk = cat(3, Hk, B);
        it=it+1;
    end    

end