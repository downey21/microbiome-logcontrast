function [filename2] = get_length_sim(p,seed,level, beta,setting, N)
n = [50, 100, 200, 500];
Suffix = {'truec','onec','noc','wrongc'};
length = zeros(p,N);

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
    length_all= get_length(result, beta, level);
    length(: , i) = length_all(2:(p+1));
    
end

filename2 = [pwd,'/result/length_sim.csv'];
csvwrite(filename2, length);

