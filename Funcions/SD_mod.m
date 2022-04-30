%{
Input:
    - H: h(x)
    - delta: valor minim que tindrem a la diagonal de la matriu de vaps
Output:

    - B: Hessiana modificada, no té cap utilitat en l'algoritme, pero a cada pas
    la guardem a la Hk
    - B1 = B^-1
    - nmod: True si no hem modificat la Hessiana
%}
function [B,Binv, nmod] = SD_mod(H,delta)
    [Q,A] = eig(H, 'vector');
    Am = max(A, delta);
    nmod = min(A == Am); % mirem si hem canviat algun eigenvalue 
    B = Q*(diag(Am))*Q'; % H simetrica => Q^-1 = Q'
    Binv = Q*diag(1./Am)*Q';
end