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
Output:
    - ak: alpha òptima
    -iWout: indicador de les condicions que compleix ak
        - 0: No WC1
        - 1: WC1
        - 2: WC2
        - 3: SWC2
        - 7: S'ha usat ELS
%}
function [a, iWout] = BLS(x, f, g, d, almin, almax, rho, c1, c2, iW, Q)
    if iW == 1 % ELS: Exact-Line Search
       a = -(g(x)'*d)/(d'*Q*d);
       iWout = 7;
    else
        a = almax;
        iWout = 0;
        [b, iWout] = WolfeC(x, a, f, g, d, c1, c2, iW);
        while a >= almin & ~b
            a = rho*a;
            [b, iWout] = WolfeC(x, a, f, g, d, c1, c2, iW);
        end
        % Per veure si falla Wolfe es pot descomentar:
        %{
        if iWout < 2
            fprintf('No solution found, iWout = %d\n', iWout);
        end
      %}  
    end
end