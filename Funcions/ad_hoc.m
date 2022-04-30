function ad_hoc
d = V(:,1);
alp = linspace(-0.2,0.2,10000);
ad_hoc = zeros(size(alp)); %inicialitzo el vector abans perque vagi més ràpid
i = 1;
for alpha = alp
    ad_hoc(i) = f(xk(:,end)+alpha*d);
    i = i+1;
end
plot(alp, ad_hoc)