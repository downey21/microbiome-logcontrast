%set path to correct directory
%check path
pwd

%import data 
load([pwd,'\IBD data analysis\IBD.mat'])
load([pwd,'\IBD data analysis\IBD_bac_name.mat'])


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
x = log(x);
x_1 = [ones(n,1) x];


%case and control id
id_control = [1:1:26];
id_case = [27:1:111];
%resampling
nrep=50;
rng(1111);
for i=1:nrep
    sample_case= sort(randsample(id_case, 56));
    sample_control = sort(randsample(id_control, 16));
    sample = [sample_control,sample_case];
    rep_set{i}.x = x_1(sample,:);
    rep_set{i}.y = y(sample);
end

% set to satisfy constraints
[n p] = size(x_1);
constr_ori = csvread([pwd,'\IBD data analysis\multi_constr.csv']);
constr=constr_ori;
[constr2, S, V] = svd(constr);
constr= constr2(: , 1:size(S,2));
Pc=constr*constr';

for i=1:nrep
x = rep_set{i}.x;
y = rep_set{i}.y;
x = x-x*Pc;
[n p] = size(x);

    for j=1:p
        norm_x(j) = 1/10*norm(x(:,j),2);
    end
    b2 = 1./norm_x;
    c =diag(b2);
    x2 =x*c;
    constr3 = c*constr;
    %constr3 = constr3';
    
    if (sum(abs(constr3))==0)
        %constr = zeros(p,1);
    else
        [constr4,S,V] = svd(constr3);
        constr3 = constr4(:,1:size(S,2));
    end
    constr3 = constr3';
    
    %useset{i}.x = x;
    useset{i}.x2 = x2;
    useset{i}.y = y;
    useset{i}.constr3 = constr3;
    useset{i}.c = c;
end

penalize =0;
result = cell(nrep);
parpool('local',4);
parfor i=1:nrep
    y = useset{i}.y;
    %x = useset{i}.x;
    x2 = useset{i}.x2;
    constr3 = useset{i}.constr3;
    c = useset{i}.c;
    Lambda=[0:0.025:0.5];
    Supp = zeros(length(Lambda),77);
    for j=1: length(Lambda)
        lambda = Lambda(j);
        [betahat,betahat_nset0]=byapg(penalize,'p',y,x2, constr3, lambda, 1, 50000, 1e-6);
        supp = betahat==0;
        supp = supp(2:78);
        Supp(j,:) = supp;
    end
    result{i} = Supp; 
end

filename = [pwd,'\IBD data analysis\Stability_multi.mat'];
save(filename, 'result');

filename = [pwd,'\IBD data analysis\sample_multi.mat'];
save(filename, 'rep_set');

delete(gcp)

Lambda=[0:0.025:0.5];
stability = zeros(length(Lambda),77);
for i=1:nrep
    stability = stability + result{i};
end
stability = stability/nrep;

filename3 = [pwd,'\IBD data analysis\Stability_multi.csv'];
csvwrite(filename3, stability);

