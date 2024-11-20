function [file1] = sim_raw_data(n, p, nsim, zeta, seed, intercept)

%set seed
rng(seed);
beta = [0.45, -0.4, 0.45, 0, -0.5, 0, 0, 0, 0, 0, -0.6, 0, 0.3, 0, 0, 0.3, zeros(1, p-16)]';% multiple constraints
beta_0 = -1;
% niu for the mean of log-normal distribution
niu = [ones(1,5)*p/2 ones(1,p-5)]';
%generate sigma matrix 
sigma = ones(p,1)*(1:p);
sigma = abs(sigma-sigma');
sigma = zeta.^sigma;

if intercept ==1  %include intercept
    i =1;
    while i <=nsim
        x = exp(mvnrnd(niu, sigma, n*2));
        x = x./(ones(p,1)*sum(x'))';
        %Z matrix
        temp_x = log(x);
        prob = temp_x*beta + beta_0;
        for j =1: (n*2)
            if prob(j)>=10
                prob1(j) =1;
            elseif prob(j) <=-10
                prob1(j) =0;
            else
                prob1(j) = exp(prob(j))/(exp(prob(j))+1);
            end
        end
        temp_y = binornd(1,prob1)'; 
        temp_x = [ones(n*2,1) temp_x];
        
        idx_case = find(temp_y==1);
        idx_control = find(temp_y ==0);
        
        if length(idx_case) >=0.4*n && length(idx_control)>=0.6*n
        sample_case= sort(randsample(idx_case, 0.4*n));
        sample_control = sort(randsample(idx_control, 0.6*n));
        %end
        sample = [sample_case;sample_control];
        data.x = temp_x(sample,:);
        data.y = temp_y(sample);
        dataset{i} = data;
        i=i+1;
        end
    end
    beta = [beta_0 ; beta];
    file1 = [pwd,'/data_n', num2str(n), '_p', num2str(p),'_seed', num2str(seed),'_sim','_zeta', num2str(zeta),'.mat'];
    save(file1 ,'beta','dataset');
else   % do not include intercept
    for i = 1:nsim
        x = exp(mvnrnd(niu, sigma, n));
        x = x./(ones(p,1)*sum(x'))';
        %Z matrix
        temp.x = log(x);
        prob = temp.x*beta;
        for j =1: n
            if prob(j)>=10
                prob1(j) =1;
            elseif prob(j) <=-10
                prob1(j) =0;
            else
                prob1(j) = exp(prob(j))/(exp(prob(j))+1);
            end
        end
        temp.y = binornd(1,prob1)'; 
        dataset{i} = temp;
    end
    file1 = [pwd,'/data/data_n', num2str(n), '_p', num2str(p), '_zeta', num2str(zeta), '_seed', num2str(seed) ,'.mat'];
    save(file1 ,'beta', 'dataset');
end


