function [betahat,betahat_nset0] = byadmm(penalize,y,x, constr2, lambda, rho, mu, default, tolerance1, tolerance2, step, step_NR)
tic
[n p] = size(x);
k = size(constr2, 1);
iter = 0;
ksi = zeros(k,1);
ksi0 = ones(k,1);
C= constr2';

if (sum(abs(constr2))==0) %no constraint
    i=0;
    beta_temp = ones(p,1)/p;
    u = zeros(p,1);
    z = zeros(p,1);
    z_temp = zeros(p,1); 
    r = 1;
    s = 1;
    while (i <= step && (r > tolerance2 | s > tolerance2))    
        %update beta_temp by Newton-Raphson method(weighted least square)
        rep=0;
        dis = 1;
        beta_temp = ones(p,1)/p;
        while (dis>default && rep <= step_NR)
            prob =x*beta_temp;
            V = diag(exp(prob)./ (exp(prob)+1).^2);
            I = x'* V * x  + rho* eye(p);
            U = -1/n * x'* (y - exp(prob)./(exp(prob)+1)) + rho* (beta_temp - z + u );
            B = I* beta_temp - U;
            temp = I \ B;
            dis = sum(abs(temp - beta_temp))/p;
            beta_temp = temp;
            rep = rep+1;
        end
        %update z and u
        for j=1:p
            z(j) = sign(beta_temp(j) + u(j))* max(0, abs(beta_temp(j) + u(j))- lambda/rho);
        end
        u = u + beta_temp - z;
        i=i+1;
        r = sum(abs(beta_temp-z));
        s = sum(abs(z-z_temp));
        z_temp =z;
    end
    if (i < step)
        disp('solved')  
    end
    betahat = beta_temp;
    betahat_nset0 =beta_temp;
    idx = find(abs(betahat) <=1e-3);
    betahat(idx) =0;
    
else
    if penalize == 0
        iter =1;
        ksi = zeros(k,1);
        ksi0 = ones(k,1);
        while (iter < step && sum(abs(ksi - ksi0)) > tolerance1)
            %update ksi
            ksi0= ksi;
    
            %setup for betahat step
            i=0;
            beta_temp = ones(p,1)/p;
            u = zeros(p,1);
            z = zeros(p,1);
            z_temp = z;
            r =1;
            s = 1;
            dis = 1;
            %betahat step
            while (i <= step && (r > tolerance2 | s > tolerance2))    
                %update beta_temp by Newton-Raphson method(weighted least square)
                rep=0;
                beta_temp = ones(p,1)/p;
                dis = 1;
                while (dis>default && rep <= step_NR)
                    prob =x*beta_temp;
                    V = diag(exp(prob)./ (exp(prob)+1).^2);
                    I = x'* V * x + mu* C*C' + rho* eye(p);
                    U = -1/n * x'* (y - exp(prob)./(exp(prob)+1)) + mu * C*C'*beta_temp + mu*C*ksi + rho* (beta_temp - z + u );
                    B = I* beta_temp - U;
                    temp = I \ B;
                    dis = norm(temp - beta_temp ,2 )/p;
                    beta_temp = temp;
                    rep = rep+1;
                end
                %update z and u
                for j=1:p
                    z(j) = sign(beta_temp(j) + u(j))* max(0, abs(beta_temp(j) + u(j))- lambda/rho);
                end
                u = u + beta_temp - z;
                i=i+1;
                r = sum(abs(beta_temp-z));
                s = sum(abs(z-z_temp));
                z_temp = z ; 
            end
            %update kasi
            ksi = ksi + C'*beta_temp;
            iter = iter+1;
        end
        if (iter < step)
            disp('solved')
        else
        %disp('un-solve')
        %disp(iter)
        end
        betahat = beta_temp;
        betahat_nset0 =beta_temp;
        idx = find(abs(betahat) <=1e-3);
        betahat(idx) =0;
        
    else  %penalize =1
        iter =1;
        ksi = zeros(k,1);
        ksi0 = ones(k,1);
        
        while (iter < step && sum(abs(ksi - ksi0)) > tolerance1)
            %update ksi
            ksi0= ksi;
    
            %setup for betahat step
            i=0;
            beta_temp = [1;ones(p-1,1)/p-1];
            %beta_temp = ones(p,1)/p;
            u = zeros(p,1);
            z = zeros(p,1);
            z_temp = z; 
            r = 1;
            s = 1;
            
            while (i <= step && (r > tolerance2 | s > tolerance2))    
                %update beta_temp by Newton-Raphson method(weighted least square)
                rep=0; dis = 1;
                beta_temp = [1;ones(p-1,1)/p-1];
                while (dis>default && rep <= step_NR)
                    prob =x*beta_temp;
                    V = diag(exp(prob)./ (exp(prob)+1).^2);
                    I = x'* V * x + mu* C*C' + rho* eye(p);
                    U = -1/n * x'* (y - exp(prob)./(exp(prob)+1)) + mu * C*C'*beta_temp + mu*C*ksi + rho* (beta_temp - z + u );
                    B = I* beta_temp - U;
                    temp = I \ B;
                    dis = norm(temp - beta_temp ,2 )/p;
                    beta_temp = temp;
                    rep = rep+1;
                end
                
                %update z and u
                z(1) = beta_temp(1) + u(1);
                for j=2:p
                    z(j) = sign(beta_temp(j) + u(j))* max(0, abs(beta_temp(j) + u(j))- lambda/rho);
                end
                u = u + beta_temp - z;
                i=i+1;
                r = sum(abs(beta_temp-z));
                s = sum(abs(z-z_temp));
                z_temp =z;
            end
            
            if i <=step
                disp('solved')
            end
            
            ksi = ksi + C'*beta_temp;
            iter = iter+1;
            %disp(ksi);
            %disp(iter);
        end
        betahat = beta_temp;
        betahat_nset0 =beta_temp;
        idx = find(abs(betahat) <=1e-3);
        betahat(idx) =0;
    end
end

toc

 
 
 
 