%{
Input:
    x1: punt inicial
    Xtr, ytr: set de training
    gL: Derivada de L
    Lte: L avaluada en el set de test
    sg_seed: Seed per les permutacions
Output:
    wo: òptim, millor punt trobat
    k: nombre d'iteracions
%}
function [wo, k] = SGM(x1, Xtr, ytr, gL,Lte, alphaSG0, betaSG, gammaSG, eSGmax, eSGbest, sg_seed)
    %xk = [x1];
    w = x1;
    wo = x1;
    p = length(Xtr);
    m = floor(gammaSG * p);
    kSGe = ceil(p / m);
    KSGmax = eSGmax * kSGe;
    e = 0;
    s = 0;
    LTEbest = inf;
    k = 0;
    kSG = floor(betaSG * KSGmax);
    alphaSG = alphaSG0 * 0.01;
    if ~isempty(sg_seed), rng(sg_seed); end 
    top = ceil(p/m-1);
    while e <= eSGmax && s < eSGbest && k < KSGmax
          P = randperm(p);
          XtrP = Xtr(:,P);
          ytrP = ytr(P);
          for i = 0:top

              XtrS = XtrP(:, (i*m+1): min((i+1)*m, p));
              ytrS = ytrP((i*m+1): min((i+1)*m, p));
              d = -gL(w, XtrS, ytrS);
              
              if k <= kSG
                  alpha_k = (1 - k/kSG)*alphaSG0 + k/kSG*alphaSG;
              else
                  alpha_k = alphaSG;
              end 
              w = w + alpha_k*d;
              %xk = [xk w];
              k = k + 1;
          end
          e = e + 1;
          LTE = Lte(w);
          if LTE < LTEbest
              LTEbest = LTE;
              wo = w;
              s = 0;
          else
              s = s + 1;
          end   
    end
    k = k + 1;
end



