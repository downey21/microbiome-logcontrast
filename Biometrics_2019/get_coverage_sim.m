function [filename2] = get_coverage_sim(p,seed,level, beta,setting, N)
n = [50, 100, 200, 500];
Suffix = {'truec','onec','noc','wrongc'};
coverage = zeros(p+1,N);

for i =1:N
    j = mod(i,4);
    if j==0
        j=4;
    end
    k = (i-j)/4+1;
    n1 = n(k);
    suffix = Suffix{j};
    filename = [pwd,'/result/res_', suffix,'_n',num2str(n1), '_p', num2str(p),'_seed',num2str(seed),setting, '.mat'];
    load(filename);
    cov_prob= get_coverage(result, beta, level);
    coverage(: , i) = cov_prob; 
    
end

filename2 = [pwd,'/result/coverage_sim.csv'];
csvwrite(filename2, coverage);

