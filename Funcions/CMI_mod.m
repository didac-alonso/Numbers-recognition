%{
Donada H(x) dona la modificació CMI de la matriu

Input:
    H: H(x)
Output:
    B: Hessiana modificada
    tau: tau usada per la modificació, tau = 0, no hem modificat
%}
function [B,tau] = CMI_mod(H)
    tau = 0;    
    laUB = norm(H, 'fro');
    k = 0;
    n = size(H,1);
    [R, err] = chol(H);
    l = 0;
    while(err > 0)
        l = l+1;
        tau = (1.01-1/2^l)*laUB;
        [R, err] = chol(H+tau*eye(n));
    end
    B = R'*R;
end