function length_CI = get_length(result, beta, level)

p = length(beta);
nsim = length(result);
length_CI = zeros(p,1);
for i = 1:nsim
    beta_u = result{i}.beta_u;
    CI_u = result{i}.CI_u; CI_l = result{i}.CI_l;
    width = (CI_u-CI_l)/2;
    width = width/norminv(0.975)*norminv(1-(1-level)/2);
    length_CI = length_CI + width*2;
end

length_CI = length_CI/nsim;
