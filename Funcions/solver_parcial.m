%{
Input:
    - isd: mètode de descens usat
        - 1 GM
        - 2 CGM
        - 3 BFGS
        - 4 NM
        - 5 MNM-SD
        - 6 MNM-CMI
    - icg:  variant de CGM (relevància si isd = 2)
        - 1 FR
        - 2 PR+
    - irc: Restart per la CGM, 
        - 0 no restart
        - 1 RC1, cada n iteracions reinci
        - 2 RC2, té a veure amb perpendicularitat
    - nu: condició de quin numero no ha de suposar RC2
    - delta: requerida pel NMN-SD, minim valor a la diagonal de la hessiana
%}
function [xk, dk, alk, iWk, betak, Hk, tauk, it] = solver_parcial(x1, f, g, h,epsG, kmax, almax, almin, rho, c1, c2, iW, isd, icg, irc, nu, delta,parcial)
    betak = [];
    Hk = [];
    tauk = [];
    Q = h(x1);
    if isd == 1
        [xk, dk, alk, iWk, it] = GM(x1, f, g, almin, almax, rho, c1, c2, iW, epsG, kmax, Q);
    elseif isd == 2 % cridem a la nova funcio
        [xk, dk, alk, betak, iWk, it] = CGM_parcial(x1, f, g, almin, almax, rho, c1, c2, iW, epsG, kmax, icg, irc, nu, Q,parcial);
    elseif isd == 3
        [xk, dk, alk, Hk, iWk, it] = BFGS(x1, f, g, almin, almax, rho, c1, c2, iW, epsG, kmax, Q);
    elseif isd == 4
        [xk, dk, alk, Hk, iWk, it] = Newton(x1, f, g, h, almin, almax, rho, c1, c2, iW, epsG,kmax);
    elseif isd == 5
        [xk, dk, alk, Hk, iWk, it] = MNM_SD(x1, f, g, h, almin, almax, rho, c1, c2, iW, epsG,kmax, delta, Q);
    elseif isd == 6
        [xk, dk, alk, Hk, iWk, tauk, it]= MNM_CMI(x1, f, g, h, almin, almax, rho, c1, c2, iW, epsG,kmax, Q);
    end
end