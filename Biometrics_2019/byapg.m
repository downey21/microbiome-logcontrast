function [betahat,betahat_nset0] = byapg(penalize,alg,y,x, constr2, lambda, mu,step, tolerance1)
%tic

[n p] = size(x);
k = size(constr2, 1);
iter = 1;
ksi = zeros(k,1);
ksi0 = ones(k,1);
C= constr2';

beta_0 = ones(p,1)/p;
y_0 = beta_0;
beta =beta_0;
beta_temp = beta_0;

r =10;
line_search = true;
gamma = 0.5;
t_k =10000;
epsilon = 1e-6;


if (sum(abs(constr2))==0)
    if penalize ==0
        y_k = y_0;
        while iter <=step
            while line_search
                prob =x*y_k;
                for j =1: n
                    if prob(j)>=10
                        prob1(j) =1;
                    elseif prob(j) <=-10
                        prob1(j) =0;
                    else
                        prob1(j) = exp(prob(j))/(exp(prob(j))+1);
                    end
                end
                
                grad_g = -1/n * x'* (y - prob1');
                u = y_k - t_k* grad_g; 
                beta_temp = sign(u) .* max(0, abs(u) - lambda* t_k);
                search = -1/n * (y'*x* beta_temp - sum(log(1+exp(x*beta_temp)))) + 1/n * (y'*x* y_k - sum(log(1+exp(x*y_k))));
                search = search +  grad_g' * (y_k - beta_temp) - t_k/2 * norm(1/t_k*(y_k - beta_temp),2)^2;
                if search >0
                    t_k = t_k* gamma;
                else
                    line_search = false;
                end
            end
            line_search = true;
            y_k = beta_temp + (iter-1)/(iter + r -1)*(beta_temp -beta);
            
            prob =x*beta_temp;
                    for j =1: n
                        if prob(j)>=10
                            prob1(j) =1;
                        elseif prob(j) <=-10
                            prob1(j) =0;
                        else
                            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
                        end
                    end
                
                    grad_g = -1/n * x'* (y - prob1');
                    u = y_k - t_k* grad_g; 
                
                
                search = -1/n * (y'*x* u - sum(log(1+exp(x*u)))) + 1/n * (y'*x* beta_temp - sum(log(1+exp(x*beta_temp))));
                search = search +  grad_g' * (beta_temp - u) - t_k/2 * norm(1/t_k*(u - beta_temp),2)^2;
                
            if abs(search)<epsilon
                betahat_nset0 = beta_temp;
                betahat = beta_temp;
                disp('solved')
                break;
            end
            
            beta = beta_temp;
            iter = iter +1;
        end
        if iter >=step
            betahat_nset0 = bycvx(penalize,y,x,constr2, lambda);
            betahat = betahat_nset0;
            disp('solved by cvx');
        end
        idx = find(abs(betahat) <=1e-3);
        betahat(idx) =0;
    else
        y_k = y_0;
        while iter <=step
            while line_search
                prob =x*y_k;
                for j =1: n
                    if prob(j)>=10
                        prob1(j) =1;
                    elseif prob(j) <=-10
                        prob1(j) =0;
                    else
                        prob1(j) = exp(prob(j))/(exp(prob(j))+1);
                    end
                end
                
                grad_g = -1/n * x'* (y - prob1');
                u = y_k - t_k* grad_g; 
                
                beta_temp(1) = u(1);
                for j=2:p
                     beta_temp(j) = sign(u(j)) .* max(0, abs(u(j)) - lambda* t_k);
                end
               
                search = -1/n * (y'*x* beta_temp - sum(log(1+exp(x*beta_temp)))) + 1/n * (y'*x* y_k - sum(log(1+exp(x*y_k))));
                search = search +  grad_g' * (y_k - beta_temp) - t_k/2 * norm(1/t_k*(y_k - beta_temp),2)^2;
                if search >0
                    t_k = t_k* gamma;
                else
                    line_search = false;
                end
            end
            line_search = true;
            y_k = beta_temp + (iter-1)/(iter + r -1)*(beta_temp -beta);
            
            prob =x*beta_temp;
                    for j =1: n
                        if prob(j)>=10
                            prob1(j) =1;
                        elseif prob(j) <=-10
                            prob1(j) =0;
                        else
                            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
                        end
                    end
                
                    grad_g = -1/n * x'* (y - prob1');
                    u = y_k - t_k* grad_g; 
                
                
                search = -1/n * (y'*x* u - sum(log(1+exp(x*u)))) + 1/n * (y'*x* beta_temp - sum(log(1+exp(x*beta_temp))));
                search = search +  grad_g' * (beta_temp - u) - t_k/2 * norm(1/t_k*(u - beta_temp),2)^2;
                
            if abs(search)<epsilon
                betahat_nset0 = beta_temp;
                betahat = beta_temp;
                break;
            end
            
            beta = beta_temp;
            iter = iter +1;
        end
         if iter >=step
            betahat_nset0 = bycvx(penalize,y,x,constr2, lambda);
            betahat = betahat_nset0;
            disp('solved by cvx');
        end
        idx = find(abs(betahat) <=1e-3);
        betahat(idx) =0;
    end
else
    if alg =='L'
     if penalize ==0
         while (iter < step && sum(abs(ksi - ksi0)) > tolerance1)
             iter_1 =1;
             y_k = y_0;
             beta =beta_0;
             beta_temp = beta_0;
             ksi0 = ksi;
             while iter_1 <=step
                while line_search
                    prob =x*y_k;
                    for j =1: n
                        if prob(j)>=10
                            prob1(j) =1;
                        elseif prob(j) <=-10
                            prob1(j) =0;
                        else
                            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
                        end
                    end
                
                    grad_g = -1/n * x'* (y - prob1') + mu *  C*C'* y_k + mu*C*ksi;
                    u = y_k - t_k* grad_g; 
                    beta_temp = sign(u) .* max(0, abs(u) - lambda* t_k);
                    search = -1/n * (y'*x* beta_temp - sum(log(1+exp(x*beta_temp))))+ mu/2 * norm(C'*beta_temp+ksi,2)^2 - mu/2 * norm(C'*y_k+ksi,2)^2 + 1/n * (y'*x* y_k - sum(log(1+exp(x*y_k))));
                    search = search +  grad_g' * (y_k - beta_temp) - t_k/2 * norm(1/t_k*(y_k - beta_temp),2)^2;
                    if search >0
                        t_k = t_k* gamma;
                    else
                        line_search = false;
                    end
                end
                line_search = true;
                y_k = beta_temp + (iter-1)/(iter + r -1)*(beta_temp -beta);
            
                prob =x*beta_temp;
                    for j =1: n
                        if prob(j)>=10
                            prob1(j) =1;
                        elseif prob(j) <=-10
                            prob1(j) =0;
                        else
                            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
                        end
                    end
                
                    grad_g = -1/n * x'* (y - prob1');
                    u = y_k - t_k* grad_g; 
                
                
                search = -1/n * (y'*x* u - sum(log(1+exp(x*u)))) + 1/n * (y'*x* beta_temp - sum(log(1+exp(x*beta_temp))));
                search = search +  grad_g' * (beta_temp - u) - t_k/2 * norm(1/t_k*(u - beta_temp),2)^2;
                
                if abs(search)<epsilon
                    betahat_nset0 = beta_temp;
                    betahat = beta_temp;
                    break;
                end
            
                beta = beta_temp;
                iter_1 = iter_1 +1;
             end
             ksi = ksi + C'*beta_temp;
             iter = iter+1;
         end
         if iter < step
             disp('solved');
         end
         betahat = beta_temp;
         betahat_nset0 =beta_temp;
         idx = find(abs(betahat) <=1e-3);
         betahat(idx) =0;
     else
         while (iter < step && sum(abs(ksi - ksi0)) > tolerance1)
             iter_1 =1;
             y_k = y_0;
             beta =beta_0;
             beta_temp = beta_0;
             ksi0 = ksi;
             while iter_1 <=step
                while line_search
                    prob =x*y_k;
                    for j =1: n
                        if prob(j)>=10
                            prob1(j) =1;
                        elseif prob(j) <=-10
                            prob1(j) =0;
                        else
                            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
                        end
                    end
                
                    grad_g = -1/n * x'* (y - prob1') + mu *  C*C'* y_k + mu*C*ksi;
                    u = y_k - t_k* grad_g; 
                    
                    beta_temp(1) = u(1);
                    for j=2:p
                         beta_temp(j) = sign(u(j)) .* max(0, abs(u(j)) - lambda* t_k);
                    end
               
                    search = -1/n * (y'*x* beta_temp - sum(log(1+exp(x*beta_temp))))+ mu/2 * norm(C'*beta_temp+ksi,2)^2 - mu/2 * norm(C'*y_k+ksi,2)^2 + 1/n * (y'*x* y_k - sum(log(1+exp(x*y_k))));
                    search = search +  grad_g' * (y_k - beta_temp) - t_k/2 * norm(1/t_k*(y_k - beta_temp),2)^2;
                    if search >0
                        t_k = t_k* gamma;
                    else
                        line_search = false;
                    end
                end
                line_search = true;
                y_k = beta_temp + (iter-1)/(iter + r -1)*(beta_temp -beta);
            
                prob =x*beta_temp;
                    for j =1: n
                        if prob(j)>=10
                            prob1(j) =1;
                        elseif prob(j) <=-10
                            prob1(j) =0;
                        else
                            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
                        end
                    end
                
                    grad_g = -1/n * x'* (y - prob1');
                    u = y_k - t_k* grad_g; 
                
                
                search = -1/n * (y'*x* u - sum(log(1+exp(x*u)))) + 1/n * (y'*x* beta_temp - sum(log(1+exp(x*beta_temp))));
                search = search +  grad_g' * (beta_temp - u) - t_k/2 * norm(1/t_k*(u - beta_temp),2)^2;
                
                if abs(search)<epsilon
                    betahat_nset0 = beta_temp;
                    betahat = beta_temp;
                    break;
                end
            
                beta = beta_temp;
                iter_1 = iter_1 +1;
             end
             ksi = ksi + C'*beta_temp;
             iter = iter+1;
         end
         if iter < step
             disp('solved');
         end
         betahat = beta_temp;
         betahat_nset0 =beta_temp;
         idx = find(abs(betahat) <=1e-3);
         betahat(idx) =0;
     end
    else
        if penalize ==0
            y_k = y_0;
            while iter < step
                while line_search
                     prob =x*y_k;
                    for j =1: n
                        if prob(j)>=10
                            prob1(j) =1;
                        elseif prob(j) <=-10
                            prob1(j) =0;
                        else
                            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
                        end
                    end
                
                    grad_g = -1/n * x'* (y - prob1');
                    u = y_k - t_k* grad_g; 
                    
                    beta_temp = sign(u) .* max(0, abs(u) - lambda* t_k);
                    %projection 
                    
                    
                    beta_proj = pinv(C*C')*C*C'*beta_temp;
                    beta_temp = beta_temp - beta_proj;
                    
                    search = -1/n * (y'*x* beta_temp - sum(log(1+exp(x*beta_temp)))) + 1/n * (y'*x* y_k - sum(log(1+exp(x*y_k))));
                    search = search +  grad_g' * (y_k - beta_temp) - t_k/2 * norm(1/t_k*(y_k - beta_temp),2)^2;
                    if search >0
                        t_k = t_k* gamma;
                    else
                        line_search = false;
                    end
                end
                line_search = true;
                y_k = beta_temp + (iter-1)/(iter + r -1)*(beta_temp -beta);     
                
                prob =x*beta_temp;
                    for j =1: n
                        if prob(j)>=10
                            prob1(j) =1;
                        elseif prob(j) <=-10
                            prob1(j) =0;
                        else
                            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
                        end
                    end
                
                    grad_g = -1/n * x'* (y - prob1');
                    u = y_k - t_k* grad_g; 
                
                
                search = -1/n * (y'*x* u - sum(log(1+exp(x*u)))) + 1/n * (y'*x* beta_temp - sum(log(1+exp(x*beta_temp))));
                search = search +  grad_g' * (beta_temp - u) - t_k/2 * norm(1/t_k*(u - beta_temp),2)^2;
                
                if abs(search)<epsilon
                    betahat_nset0 = beta_temp;
                    betahat = beta_temp;
                    %
                    disp('solved')
                    break;
                end
            
                beta = beta_temp;
                iter = iter +1;
            end
        if iter >=step
            betahat_nset0 = bycvx(penalize,y,x,constr2, lambda);
            betahat = betahat_nset0;
            disp('solved by cvx')
        end
        idx = find(abs(betahat) <=1e-3);
        betahat(idx) =0;
        else
            y_k = y_0;
            while iter < step
                while line_search
                     prob =x*y_k;
                    for j =1: n
                        if prob(j)>=10
                            prob1(j) =1;
                        elseif prob(j) <=-10
                            prob1(j) =0;
                        else
                            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
                        end
                    end
                
                    grad_g = -1/n * x'* (y - prob1');
                    u = y_k - t_k* grad_g; 
                    
                    beta_temp(1) = u(1);
                    for j=2:p
                         beta_temp(j) = sign(u(j)) .* max(0, abs(u(j)) - lambda* t_k);
                    end
                    
                    %projection
                    beta_proj = pinv(C*C')*C*C'*beta_temp;
                    beta_temp = beta_temp - beta_proj;
                    
                    search = -1/n * (y'*x* beta_temp - sum(log(1+exp(x*beta_temp)))) + 1/n * (y'*x* y_k - sum(log(1+exp(x*y_k))));
                    search = search +  grad_g' * (y_k - beta_temp) - t_k/2 * norm(1/t_k*(y_k - beta_temp),2)^2;
                    if search >0
                        t_k = t_k* gamma;
                    else
                        line_search = false;
                    end
                end
                line_search = true;
                y_k = beta_temp + (iter-1)/(iter + r -1)*(beta_temp - beta);
                
                
                 prob =x*beta_temp;
                    for j =1: n
                        if prob(j)>=10
                            prob1(j) =1;
                        elseif prob(j) <=-10
                            prob1(j) =0;
                        else
                            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
                        end
                    end
                
                    grad_g = -1/n * x'* (y - prob1');
                    u = y_k - t_k* grad_g; 
                
                
                search = -1/n * (y'*x* u - sum(log(1+exp(x*u)))) + 1/n * (y'*x* beta_temp - sum(log(1+exp(x*beta_temp))));
                search = search +  grad_g' * (beta_temp - u) - t_k/2 * norm(1/t_k*(u - beta_temp),2)^2;
                
                if abs(search)<epsilon
                    betahat_nset0 = beta_temp;
                    betahat = beta_temp;
                    break;
                end
            
                beta = beta_temp;
                iter = iter +1;
            end
            if iter<step
                %disp(iter);
                %disp('solved');
            end
            if iter >=step
                betahat_nset0 = bycvx(penalize,y,x,constr2, lambda);
                betahat = betahat_nset0;
                disp('solved by cvx')
            end
        idx = find(abs(betahat) <=1e-3);
        betahat(idx) =0;
        end
    end
end

idx = find(abs(betahat) <=1e-3);
betahat(idx) =0;

%toc









