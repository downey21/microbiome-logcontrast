%prepare files for making plots

%set path to correct directory
%check path
pwd

%import data
load([pwd,'\IBD data analysis_sensitivity\IBD.mat'])
y = [zeros(26,1); ones(85,1)];
w = IBD1>0;
ind_d = sum(w);


%pick main bacteria
pick = find(ind_d>= 111*0.2);
x_or = IBD1(: , pick);

x_min=min( x_or(x_or>0));
% eliminate 0
x_or = x_or + 0.5*x_min-min(x_or, 0.5*x_min);

%convert to compositional data
[n p] = size(x_or);
x = x_or./(ones(p,1)*sum(x_or'))';
x_porp = x;
x = log(x);
x_1 = [ones(n,1) x];

%import results
load([pwd,'\IBD data analysis\results\IBD_result_vs.mat']);

%score
score = x_1 * res.beta_u;
score = [score,y];
file = [pwd,'\IBD data analysis\results\score_prob.csv'];
csvwrite(file, score);

%score_five
variable = ((res.CI_u < 0)|(res.CI_l > 0))*1;
variable_id = find(variable==1);
beta_u_5 = res.beta_u([variable_id]);
x_five = x_1(:,variable_id);

score_five = x_five * beta_u_5;
score_five = [score_five,y];
file = [pwd,'\IBD data analysis\results\score_prob_5.csv'];
csvwrite(file, score_five);

%individual plot
variable = ((res.CI_u < 0)|(res.CI_l > 0))*1;
variable_id = find(variable==1);
variable_id = variable_id(variable_id>1)-1;
pick_id = pick(variable_id);


x_select = x_porp(:,pick_id);
file = [pwd,'\IBD data analysis\results\x_select.csv.csv'];
csvwrite(file, x_select);


%rank test
file = [pwd,'\IBD data analysis\log_proportion.csv.csv'];
csvwrite(file, x);






