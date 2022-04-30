%{
Utilitza el mètode quasi Newton BFGS:
    On estimem la hessiana per H = (I-p*s*y')*H*(I-p*y*s')+p*(s*s')
    Amb
        y = dfx-dfx1
        s = x-xk(:,end)
        p = 1/(y'*s)
Output:
    xk,dk,ak,iWk,it: com sempre
    Hk: vector de matrius H en cada iteració
%}
function [xk, dk, ak, Hk, iWk, it] = BFGS(x, f, g, almax, c1, c2, epsG, kmax, ialmax, maxiter, epsal)
    I = eye(size(x,1));
    H = I;
    Hk = [H]; it = 1; xk = [x]; dk = []; ak = []; iWk = [];
    dfx = g(x);
    while norm(dfx) > epsG & it <= kmax
        d = -H*dfx;
        [a,iout] = uo_BLSNW32(f,g,x,d,almax,c1,c2,maxiter,epsal);
        x = x+a*d;
        
        dfx1 = dfx;
        dfx = g(x);
        s = x-xk(:,end); y = dfx-dfx1; 
        p = 1/(y'*s);
        H = (I-p*s*y')*H*(I-p*y*s')+p*(s*s');

        if ialmax == 1 && it > 1
            almax = a*(dfx1'*dk(:,end))/(dfx'*d);
        elseif ialmax == 2
            almax = 2*(f(x)-f(xk(:,end)))/(dfx'*d);
        end
        
        
        xk = [xk x];
        dk = [dk d];
        ak = [ak a];
        iWk = [iWk iout];
        Hk = cat(3, Hk, H); % concatenacio 3D
        it = it + 1;
    end
end