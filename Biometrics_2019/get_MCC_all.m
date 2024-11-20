function [filename2] = get_MCC_all(p,seed,level, beta,setting, N)
n = [50,100, 200, 500];
Suffix = {'truec','onec','noc','wrongc'};
prediction = zeros(length(n),1+4*2);
prediction(1:(length(n)),1) = n;
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
    [TPR, FPR] = get_MCC(result, beta, level);
    prediction(k,((j-1)*2+2):((j-1)*2+3)) =  [TPR, FPR];
end

filename2 = [pwd,'/result/prediction_sim',num2str(p),'.csv'];
csvwrite(filename2, prediction);

