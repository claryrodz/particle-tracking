function trunc_dist = levy_distro_exptrunc(alpha,bet,gam,trunc,N)
alpha=alpha+1;
trunc_dist = zeros(N,1);
for i = 1:N
    g=0;
    while g==0
        z = stblrnd(alpha-1,bet,gam,0,1,1);
        uniform = rand;
        expcomparison = exp(-abs(z)/trunc);
        if expcomparison>uniform
            g = 1;
        end
    end
    trunc_dist(i) = z;
end

