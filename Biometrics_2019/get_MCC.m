function [TPR, FPR] = get_MCC(result, beta, level)

p = length(beta);
nsim = length(result);
P = zeros(nsim,1);
FP = zeros(nsim,1);
TP = zeros(nsim,1);
FN = zeros(nsim,1);
TN = zeros(nsim,1);
RP = beta~=0;
RP= RP(2:p);

for i = 1:nsim
    beta_u = result{i}.beta_u;
    CI_u = result{i}.CI_u; CI_l = result{i}.CI_l;
    width = (CI_u-CI_l)/2;
    width = width/norminv(0.975)*norminv(1-(1-level)/2);
    result{i}.CI_u = beta_u + width;
    result{i}.CI_l = beta_u - width;
    PP = (result{i}.CI_u<0)|(result{i}.CI_l>0);
    PP = PP(2:p);
    FP(i) = sum(PP&(~RP));
    TP(i) = sum(PP&RP);
    FN(i) = sum((~PP)&RP);
    TN(i) = sum((~PP)&(~RP));
    P(i) = sum(PP);
end

TPR = mean(TP./(TP+FN));
FPR = mean(FP./(FP+TN));



