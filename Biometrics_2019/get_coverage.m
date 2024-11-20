function cov_prob = get_coverage(result, beta, level)


p = length(beta);
nsim = length(result);
cov_prob = zeros(p,1);
for i = 1:nsim
    beta_u = result{i}.beta_u;
    CI_u = result{i}.CI_u; CI_l = result{i}.CI_l;
    width = (CI_u-CI_l)/2;
    width = width/norminv(0.975)*norminv(1-(1-level)/2);
    result{i}.CI_u = beta_u + width;
    result{i}.CI_l = beta_u - width;    
    cov_prob = cov_prob + ((result{i}.CI_u > beta)&(result{i}.CI_l < beta))*1;
end

cov_prob = cov_prob/nsim;



