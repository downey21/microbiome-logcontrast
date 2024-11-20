%set to correct working dirctory
%check path
pwd

%import data 
load([pwd,'\IBD data analysis_sensitivity\IBD.mat'])
load([pwd,'\IBD data analysis_sensitivity\IBD_bac_name.mat'])

y = [zeros(26,1); ones(85,1)];
w = IBD1>0;
ind_d = sum(w);


%pick main bacteria
pick = find(ind_d>= 111*0.2);
x_or = IBD1(: , pick);

x_min=min( x_or(x_or>0));
% eliminate 0
x_or = x_or + 0.1*x_min-min(x_or, 0.1*x_min);
%
% sensitivity analysis of eliminating zeros  
%


%convert to compositional data
[n p] = size(x_or);
x = x_or./(ones(p,1)*sum(x_or'))';
x = log(x);
x = [ones(n,1) x];

train_set_y = [y(1:16);y(27:82)];
test_set_y = [y(17:24);y(83:110)];

train_set_x = [x(1:16,:);x(27:82,:)];
test_set_x = [x(17:24,:);x(83:110,:)];


%% onec variable selection
[n p] = size(x);
constr =  [0;ones(p-1,1)];
[constr2, S, V] = svd(constr);
constr= constr2(: , 1:size(S,2));
Pc=constr*constr';
x = x - x*Pc;

for j=1:p
    norm_x(j) = 1/10*norm(x(:,j),2);
end
b2 = 1./norm_x;
c =diag(b2);
x_one =x*c;
constr3 = constr'*c;
[constr2, S, V] = svd(constr3');
constr3= constr2(: , 1:size(S,2));
constr3= constr3';

mu=1; level =0.95; penalize =0; length1 = 15; length2 =25;

d=1.5;
[res.beta_n,res.lambda_best, res.EBIC] = biased_estimate_BIC(penalize,y,x_one,constr3, mu, length1, length2);
[res.beta_u,res.CI_l, res.CI_u, res.CI_M] = debiased_cvx(y,x_one, res.beta_n, constr3, res.lambda_best, level,c,d);
res.beta_n = c*res.beta_n;
sum(res.beta_n(2:78))+sum(res.beta_u(2:78))


width = (res.CI_u-res.CI_l)/2;
width = width/norminv(0.975)*norminv(1-(1-level)/2);
% find the correct bacteria;
variable = ((res.CI_u < 0)|(res.CI_l > 0))*1;
variable_id = find(variable==1);
variable_id = variable_id(variable_id>1)-1;
pick_id = pick(variable_id);
bacteria = IBD_bac_name(1,pick_id);

CI = zeros(78*2,2);
CI(1:78,1) = res.beta_u;
CI(79:(78*2),1) = res.beta_n;

CI(1:78,2) = (res.CI_u-res.CI_l)/2;

%%