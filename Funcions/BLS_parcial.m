%{
Input:
    - x: punt actual
    - d: direcció de descens
    - c1, c2: per Wolfe
    - p: rho, com fem descendir alpha
    - a: alpha, donada per interval amax, amin
    - iW: metode Wolfe a usar
        - 0: ELS
        - 1: WC
        - 2: SWC
    - parcial:
        - 0: no inclou modificació
        - 1: inclou modificacio
Output:
    - ak: alpha òptima
    -iWout: indicador de les condicions que compleix ak
        - 0: No WC1
        - 1: WC1
        - 2: WC2
        - 3: SWC2
        - 7: S'ha usat ELS
%}
function [a, iWout] = BLS_parcial(x, f, df, d, amin, amax, p, c1, c2, iW, Q, parcial)
    if iW == 0 % ELS: Exact-Line Search
       a = -(df(x)'*d)/(d'*Q*d);
       iWout = 7;
    else
        a = amax;
        iWout = 0;
        [b, iWout] = WolfeC(x, a, f, df, d, c1, c2, iW);
        if b & parcial == 1 % Sol cal comprovar nova condició si es compleix WC
            if df(x+a*d)'*d >= 0 % Nova condicio negada
                b = 0;
            end
        end
        while a >= amin & ~b
            a = p*a;
            [b, iWout] = WolfeC(x, a, f, df, d, c1, c2, iW);
            if b & parcial == 1 % Sol cal comprovar nova condició si es compleix WC
               if df(x+a*d)'*d >= 0 % Nova condició negada
                  b = 0;
                end 
            end
        end
        % Per veure si falla Wolfe es pot descomentar:
        %{
        if iWout < 2
            fprintf('No solution found, iWout = %d\n', iWout);
        end
      %}  
    end
end