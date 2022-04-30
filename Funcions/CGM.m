%{
Usa el mètode del gradient conjugat per buscar punts crítics:
Aquest mètode consisteix en fer que: d = -dfx + B*dk1
    Pot usar el FR on: B = (dfx'*dfx)/(norm(dfx1))^2
    O el PR+ on B = max(0, (dfx'*(dfx-dfx1))/(norm(dfx1))^2)
        PR+ pot tenir problemes de convergència si no es compleix la
        SDC, suficient descent condition
            
Input:
    - x: Punt inicial
    - f, df: funció i derivada
    - amin, amax: rang per les alphas, s'usa en el BLS
    - c1, c2: constants per les Wolfe Conditions
    - iW: Condicions Wolfe a satisfer:
    - tol: precisió buscada
    - itmax: màxim d'iteracions que pot dur a terme l'algoritme
    - icg:  variant de CGM icg = 1 FR, icg = 2 PR+
    - irc: Restart per la CGM, irc = no, irc = 1 RC1, irc = 2 RC2
    - nu: nombre de iteracions entre restarts en cas de que irc = 1

Output:
    - x: Punt que fa 0 df.
    - it: nombre d'iteracions usades
%}
function [xk, dk, ak, Bk, iWk, it] = CGM(x, f, df, amin, amax, p, c1, c2, iW, tol, itmax, icg, irc, nu, Q)
    it = 0;
    xk = [x]; Bk = [0]; ak = []; iWk = [];
    dfx = df(x); dk = [-dfx]; d = -dfx; 
    while norm(dfx) > tol & it < itmax
        [a, iWout] = BLS(x,f,df,d,amin,amax, p, c1, c2, iW, Q);
        x = x + a*d;
        dfx1 = dfx; % df k-1
        dfx = df(x); % df k
        if irc == 1 & mod(it,length(x)) == 0 %RC1: Si it = nk, on n es dm(x) 
            B = 0;
        elseif irc == 2 & abs(dfx'*dfx1)/(norm(dfx))^2 >= nu %RC2
            B = 0;
        elseif icg == 1 % FR
            B = (dfx'*dfx)/(norm(dfx1))^2;
        elseif icg == 2 % PR+
            B = max(0, (dfx'*(dfx-dfx1))/(norm(dfx1))^2);   
        end
        d = -dfx + B*dk(:,end);

        iWk = [iWk iWout];
        xk = [xk x];
        ak = [ak a];
        Bk = [Bk B];
        dk = [dk d];
        it = it + 1;
    end
    it = it+1;
end