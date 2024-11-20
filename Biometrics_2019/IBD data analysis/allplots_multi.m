%prepare files for making plots using multiple constraints

%set path to correct directory
%check path
pwd

%import data
load([pwd,'\IBD data analysis\IBD.mat'])
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
load([pwd,'\IBD data analysis\results\IBD_result_vs_multi.mat']);

%score
score = x_1 * res.beta_u;
score = [score,y];
file = [pwd,'\IBD data analysis\results\score_prob_multi.csv'];
csvwrite(file, score);

%score_six
variable = ((res.CI_u < 0)|(res.CI_l > 0))*1;
variable_id = find(variable==1);
beta_u_6 = res.beta_u([variable_id]);
x_six = x_1(:,variable_id);

score_six = x_six * beta_u_6;
score_six = [score_six,y];
file = [pwd,'\IBD data analysis\results\score_prob_6_multi.csv'];
csvwrite(file, score_six);

%individual plot
%renormalizing data
constr_ori = csvread([pwd,'\IBD data analysis\multi_constr.csv']);
[~,r]=size(constr_ori);
x_porp_renormal = x_porp;
for i=1:r
    con = constr_ori(2:78,i);
    con = logical(con);    
    x_temp = x_porp(:,con);
    [~,p_temp] = size(x_temp);
    x_temp = x_temp./(ones(p_temp,1)*sum(x_temp'))';
    x_porp_renormal(:,con) = x_temp;
end

variable = ((res.CI_u < 0)|(res.CI_l > 0))*1;
variable_id = find(variable==1);
variable_id = variable_id(variable_id>1)-1;
pick_id = pick(variable_id);

x_select = x_porp_renormal(:,pick_id);
file = [pwd,'\IBD data analysis\results\x_select_multi.csv'];
csvwrite(file, x_select);






