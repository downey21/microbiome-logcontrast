function[betahat_ebic, lambda_best_ebic, EBIC] = biased_estimate_BIC(penalize,y,x,constr2, mu,length1,length2)

%calcuate the biased estimate by admm and cvx. 
%tic
[n p] = size(x);
EBIC = zeros(5,p);
gamma = 1 - 1/(2*log(p)/log(n));

%gamma =0;
fin=0;
count = [0,1];
lambda =[0];
lambda0 = 1/8;  %first try
lambda_l= 0; 
lambda_temp = lambda0;
beacon =1;
bound = [ones(1,p+1);zeros(1,p+1)]; % lower and upper bound
search_step =1;
num_count = 1;
num_value = 1;

while (fin ==0 && search_step <=500) 
   lambda_try = lambda_temp;
   
   %[betahat, betahat_nset0] = byadmm(penalize,y,x, constr2, lambda_try, 1, 1, 1e-6, 1e-5, 1e-6, 20000, 1000);
   %[betahat,betahat_nset0] = byapg(penalize,'p',y,x, constr2, lambda_try, mu, 50000, 1e-6);
   if beacon <=3
        betahat_nset0 = bycvx(penalize,y,x,constr2, lambda_try);
        betahat = betahat_nset0;
        idx = find(abs(betahat) <=1e-3);
        betahat(idx) =0;
   else
       [betahat,betahat_nset0] = byapg(penalize,'p',y,x, constr2, lambda_try, mu, 50000, 1e-6);
   end
   
   niu = sum(betahat~=0);
       
   if lambda_try <=bound(1,niu+1)
       bound(1,niu+1) = lambda_try;
   end
   if lambda_try >= bound(2,niu+1)
       bound(2,niu+1) = lambda_try;
   end
   
   if beacon<=max(count)
       idx = min(find(count>=beacon));
       lambda_l = bound(2, count(idx)+1);
   else
       lambda_l=0;
   end
   
   if niu >= beacon
       nlik = -2*(y'*x*betahat_nset0 - sum(log(1+exp(x*betahat_nset0))));
       ebic = nlik + niu * log(n) + 2* niu * gamma * log(p);
       bic = nlik+ niu*log(n);
       if any(count==niu)==0
           count = [count, niu];           
           EBIC(1,p+1-niu) = ebic;
           EBIC(2,p+1-niu) = bic;
           EBIC(3,p+1-niu) = nlik;
           EBIC(4,p+1-niu) = lambda_try;
           EBIC(5,p+1-niu) = niu;
       else
           if EBIC(1,p+1-niu)>= ebic
               EBIC(1,p+1-niu) = ebic;
               EBIC(4,p+1-niu) = lambda_try;
           end
       end
      
       %reset beacon
       count = sort(count);
       search =beacon;
       stop=0;
       while stop ==0
           if any(count == search)==1
               search = search +1;
           else
               stop =1;
           end
       end
       beacon =search;
       %new lambda_temp
       %if beacon >=2
           lambda_temp = (lambda_try + bound(1,beacon))/2;
       %else
           %lambda_temp = lambda_try*2;
       %end    
   else
       lambda_temp = (lambda_try + lambda_l)/2;
   end
   
   [~,num] =size(count);
   
   if num == num_value
       num_count = num_count +1;
   else
       num_value = num;
       num_count =0;
   end
   
   if num_count >= (50 - num)
       count = sort([count, beacon]);
       num_count =0;
   end
  
  
   disp([count,search_step])
   lo=0;
   idx1 = find(EBIC(1,1:p)>0);
   if beacon >=length1
       if (EBIC(1, idx1(1)) > EBIC(1, idx1(2)) && EBIC(1, idx1(2)) > EBIC(1, idx1(3)))
       lo=1;
       end
   end
   
   if (beacon >= length2 || lo==1)
       fin =1;
   end
   search_step = search_step +1;
end

idx1 = find(EBIC(1,1:p)>0);
[~, idx2] = min(EBIC(1,idx1));
lambda_best_ebic = EBIC(4, idx1(idx2));

idx = find(EBIC(2,1:p)>0);
[~, idx3] = min(EBIC(2,idx));
lambda_best_bic = EBIC(4, idx(idx3));

[~,betahat_ebic] = byapg(penalize,'P',y,x, constr2, lambda_best_ebic, 1, 50000, 1e-6);



%toc


