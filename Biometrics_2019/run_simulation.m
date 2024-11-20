% cd ... set to correct folder
%%parameter setup
nsim =100;
%rho =1; mu =1; default =1e-6; tolerance1 = 1e-4; tolerance2 = 1e-6;
%step =1000; step_NR=100;
level =0.95;

%% case1
intercept =1; penalize =0;
p=50; zeta =0.2;
constr = get_constr(p,1,intercept);
constr0 = zeros(p+1,1);
constr1 = [0;ones(p,1)];
constr2 = get_wrong_constr(p);
setting = '_sim_zeta0.2';
seed =322;

n=200;
sim_raw_data(n, p, nsim, zeta, seed, intercept);
parallel(penalize,n,p,nsim, seed, setting, constr,1,level,'truec');
parallel(penalize,n,p,nsim, seed, setting, constr1,1,level,'onec');
parallel(penalize,n,p,nsim, seed, setting, constr0,1,level,'noc');
parallel(penalize,n,p,nsim, seed, setting, constr2,1,level,'wrongc');
n=100;
sim_raw_data(n, p, nsim, zeta, seed, intercept);
parallel(penalize,n,p,nsim, seed, setting, constr,1,level,'truec');
parallel(penalize,n,p,nsim, seed, setting, constr1,1,level,'onec');
parallel(penalize,n,p,nsim, seed, setting, constr0,1,level,'noc');
parallel(penalize,n,p,nsim, seed, setting, constr2,1,level,'wrongc');
n=50;
sim_raw_data(n, p, nsim, zeta, seed, intercept);
parallel(penalize,n,p,nsim, seed, setting, constr,1,level,'truec');
parallel(penalize,n,p,nsim, seed, setting, constr1,1,level,'onec');
parallel(penalize,n,p,nsim, seed, setting, constr0,1,level,'noc');
parallel(penalize,n,p,nsim, seed, setting, constr2,1,level,'wrongc');
n=500;
sim_raw_data(n, p, nsim, zeta, seed, intercept);
parallel(penalize,n,p,nsim, seed, setting, constr,1,level,'truec');
parallel(penalize,n,p,nsim, seed, setting, constr1,1,level,'onec');
parallel(penalize,n,p,nsim, seed, setting, constr0,1,level,'noc');
parallel(penalize,n,p,nsim, seed, setting, constr2,1,level,'wrongc');

%% coverage probability and length of CI
N = 16;
beta = [-1,0.45, -0.4, 0.45, 0, -0.5, 0, 0, 0, 0, 0, -0.6, 0, 0.3, 0, 0, 0.3, zeros(1, p-16)]';
get_coverage_sim(p,seed,level, beta,setting, N);
get_length_sim(p,seed,level, beta,setting, N);

%% prediction evaluation
get_MCC_all(p,seed,level, beta,setting, N);

%% case 2
intercept =1; penalize =0;
p=100; zeta =0.2;
constr = get_constr(p,1,intercept);
constr0 = zeros(p+1,1);
constr1 = [0;ones(p,1)];
constr2 = get_wrong_constr(p);
setting = '_sim_zeta0.2';
seed =322;

n=200;
sim_raw_data(n, p, nsim, zeta, seed, intercept);
parallel(penalize,n,p,nsim, seed, setting, constr,1,level,'truec');
parallel(penalize,n,p,nsim, seed, setting, constr1,1,level,'onec');
parallel(penalize,n,p,nsim, seed, setting, constr0,1,level,'noc');
parallel(penalize,n,p,nsim, seed, setting, constr2,1,level,'wrongc');
n=100;
sim_raw_data(n, p, nsim, zeta, seed, intercept);
parallel(penalize,n,p,nsim, seed, setting, constr,1,level,'truec');
parallel(penalize,n,p,nsim, seed, setting, constr1,1,level,'onec');
parallel(penalize,n,p,nsim, seed, setting, constr0,1,level,'noc');
parallel(penalize,n,p,nsim, seed, setting, constr2,1,level,'wrongc');
n=500;
sim_raw_data(n, p, nsim, zeta, seed, intercept);
parallel(penalize,n,p,nsim, seed, setting, constr,1,level,'truec');
parallel(penalize,n,p,nsim, seed, setting, constr1,1,level,'onec');
parallel(penalize,n,p,nsim, seed, setting, constr0,1,level,'noc');
parallel(penalize,n,p,nsim, seed, setting, constr2,1,level,'wrongc');
n=50;
sim_raw_data(n, p, nsim, zeta, seed, intercept);
parallel(penalize,n,p,nsim, seed, setting, constr,1,level,'truec');
parallel(penalize,n,p,nsim, seed, setting, constr1,1,level,'onec');
parallel(penalize,n,p,nsim, seed, setting, constr0,1,level,'noc');
parallel(penalize,n,p,nsim, seed, setting, constr2,1,level,'wrongc');
%% coverage probability and length of CI
N = 16;
beta = [-1,0.45, -0.4, 0.45, 0, -0.5, 0, 0, 0, 0, 0, -0.6, 0, 0.3, 0, 0, 0.3, zeros(1, p-16)]';
get_coverage_sim(p,seed,level, beta,setting, N);
get_length_sim(p,seed,level, beta,setting, N);
%% prediction evaluation
get_MCC_all(p,seed,level, beta,setting, N);

