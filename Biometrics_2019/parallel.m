function [filename2] = parallel(penalize,n,p,nsim, seed, setting, constr,mu,level,suffix)

filename =  [pwd,'/simulation data/data_n', num2str(n), '_p', num2str(p), '_seed', num2str(seed), setting ,'.mat'];
filename2 = [pwd,'/simulation result/res_', suffix,'_n',num2str(n), '_p', num2str(p),'_seed',num2str(seed),setting, '.mat'];
load(filename);   

if (sum(abs(constr))==0)
        constr5 = constr;
else
        [constr4,S,V] = svd(constr);
        constr5 = constr4(:,1:size(S,2));
end
Pc = constr5*constr5';

 length1 = 10;
 length2 = 15;

 d = 0.01;
 
for i=1:nsim
x = dataset{i}.x;
y = dataset{i}.y;
x = x-x*Pc;
[n p] = size(x);



    for j=1:p
        norm_x(j) = 1/10*norm(x(:,j),2);
    end
    b2 = 1./norm_x;
    c =diag(b2);
    x2 =x*c;
    constr3 = c*constr5;
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

result = cell(nsim);

parpool('local',4); %change to the number of cores of specific machine. 
parfor i=1:nsim
    y = useset{i}.y;
    %x = useset{i}.x;
    x2 = useset{i}.x2;
    constr3 = useset{i}.constr3;
    c = useset{i}.c;
    [result{i}.beta_n,result{i}.lambda_best, result{i}.EBIC] = biased_estimate_BIC(penalize,y,x2,constr3, 1,length1,length2);
    [result{i}.beta_u,result{i}.CI_l,result{i}.CI_u, result{i}.CI_M] = debiased_cvx(y,x2, result{i}.beta_n, constr3, result{i}.lambda_best, level,c,d);
    %[res.beta_u res.CI_l res.CI_u res.CI_M] = debiased_admm(y,x, res.beta_n, constr2,  res.lambda_best, level,tolerance2, step);
    result{i}.beta_n = c*result{i}.beta_n; 
    result{i}.check = sum(result{i}.beta_u(2:p)) + sum(result{i}.beta_n(2:p));
end


save(filename2, 'beta', 'result');
delete(gcp)