%{
Busca zeros a la derivada d'una funció.

Input:
    - x: Punt inicial
    - f, g: funció i derivada
    - almin, almax: rang per les alphas, s'usa en el BLS
    - c1, c2: constants per les Wolfe Conditions
    - iW: Condicions Wolfe a satisfer, 
        - iW= 0 ELS
        - iW = 1 WC
        - iW = 2 SWC
    - tol: precisió buscada
    - itmax: màxim d'iteracions que pot dur a terme l'algoritme
    - Q: s'usa en cas de funció quadràtica que usi el ELS

Output:
    - xk: 
    - it: nombre d'iteracions usades
%}
function [xk, dk, ak, iWk, it] = GM(x, f, g, almin, almax, c1, c2, epsG, itmax,ialmax, maxiter,epsal, ils, Q)
    it = 1; xk = [x]; dfx = g(x); rho = 0.05;
    dk = [-dfx]; d = -dfx; ak = []; iWk = [];
    while norm(dfx) > epsG && it <= itmax
        if ils == 3 % BLS from N&W
            [a,iout] = uo_BLSNW32(f,g,x,d,almax,c1,c2,maxiter,epsal);            
        elseif ils < 3 % BLS o ELS
            [a, iout] = BLS(x, f, df, d, almin, almax, rho, c1, c2, iW, Q);
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