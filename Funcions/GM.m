%{
Busca zeros a la derivada d'una funció.

Input:
    - xk: Punt inicial
    - f, df: funció i derivada
    - amin, amax: rang per les alphas, s'usa en el BLS
    - c1, c2: constants per les Wolfe Conditions
    - iW: Condicions Wolfe a satisfer, 
        - iW= 0 ELS
        - iW = 1 WC
        - iW = 2 SWC
    - tol: precisió buscada
    - itmax: màxim d'iteracions que pot dur a terme l'algoritme
    - Q: s'usa en cas de funció quadràtica que usi el ELS

Output:
    - xk,dk,ak, iWk com sempre
    - it: nombre d'iteracions usades
%}
function [xk, dk, ak, iWk, it] = GM(x, f, g, almax, c1, c2, epsG, itmax,ialmax, maxiter,epsal, ils, Q)
    it = 1;
    xk = [x];
    dfx = g(x);
    dk = [-dfx]; d = -dfx; ak = []; iWk = [];
    while norm(dfx) > epsG && it <= itmax
        if ils == 3
            [a,iout] = uo_BLSNW32(f,g,x,d,almax,c1,c2,maxiter,epsal);            
        elseif ils < 3
            [a, iout] = BLS(x, f, df, d, amin, amax, p, c1, c2, iW, Q);
        end
        x = x + a*d;
        dfx = g(x);
        d = - dfx;
        if ialmax == 1
            almax = a*(g(xk(:,end))'*dk(:,end))/(dfx'*d);
        elseif ialmax == 2
            almax = 2*(f(x)-f(xk(:,end)))/(dfx'*d);
        end
        xk = [xk x];
        ak = [ak a];
        dk = [dk d];
        iWk = [iWk iout];
        it = it + 1;
    end
end