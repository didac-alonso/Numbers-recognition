%{
Busca zeros a la derivada d'una funció.

Input:
    - xk: Punt inicial
    - f, df, d2f: funció, derivada i segona derivada
    - amin, amax: rang per les alphas, s'usa en el BLS
    - c1, c2: constants per les Wolfe Conditions
    - iW: Condicions Wolfe a satisfer:
    - tol: precisió buscada
    - itmax: màxim d'iteracions que pot dur a terme l'algoritme

Output:
    - xk: Punt que fa 0 df.
    - it: nombre d'iteracions usades
%}
function [xk, dk, alk, Hk, iWk, it]=Newton(x, f, df, d2f, amin, amax, p, c1, c2, iW, tol,itmax)
    it=1;
    xk = [x];
    dfk = df(x);
    Hk = [];
    a = 1;
    dk = [];
    while norm(dfk)>tol & it<=itmax
        H = d2f(x);
        d = -H\dfk; % Equivalent a d2f(xk)^-1*dfk, però és més eficient així
        x=x+a*d;
        dfk = df(x);
        
        xk = [xk x];
        dk = [dk d];
        Hk = cat(3, Hk, H);
        it=it+1;
        
    end    
    alk = 1+zeros(it);
    iWk = alk + 3;
end