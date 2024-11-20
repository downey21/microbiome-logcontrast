function [beta_u, CI_l, CI_u, M, fail] = debiased_cvx(y,x, betahat, constr, lambda, level,c,d)

[n p] = size(x);

fail = 0;

prob =x*betahat;
V = diag(exp(prob)./ (exp(prob)+1).^2);
Sigma = x'* V *x /n;
Q = eye(p) - constr'*constr;

[U_sig,S_sig,V_sig] = svd(Sigma);
S_sig = diag(S_sig);
eps = S_sig(size(S_sig,1)-1)/1000;
Sig_eps = Sigma + eps*eye(p);

%all simulation choice;
% try 
gamma = d*lambda;


M = zeros(p,p);
b1 = ones(p,1)*gamma;



for i = 1:p
    %disp(i)
    b2 = Q(:,i);
    low = b2 - b1;
    up = b2 + b1;
   
    %find m with initial gamma
    cvx_begin quiet
    cvx_precision low
    variable m(p)
    minimize(m' * Sig_eps * m)
        subject to
        low <= Sigma*m <= up
	cvx_end
    M(i, :) = m;
    
    % if not solved, need to update gamma
    if(~(strcmp(cvx_status,'Solved')|| strcmp(cvx_status,'Inaccurate/Solved')))
        disp('first_failed')
        cvx_begin quiet
        cvx_precision low
            variables gamma2 m(p)
            minimize(gamma2)
            subject to
            b2 - ones(p,1)*gamma2 <= Sigma*m <= b2 + ones(p,1)*gamma2
        cvx_end
        M(i,:) = m; 
        low2 = b2 - ones(p,1)*max(1.1*gamma2,1.2*gamma);
        up2 = b2 + ones(p,1)*max(1.1*gamma2,1.2*gamma);

        cvx_begin quiet
        cvx_precision low
            variables m(p)
            minimize(m'*Sig_eps*m)
            subject to
            low2 <= Sigma*m <= up2
        cvx_end
        if( strcmp(cvx_status,'Solved')|| strcmp(cvx_status,'Inaccurate/Solved') )
            M(i,:) = m;
            disp('second_solved')
        else
            fail =fail+1;
        end
    end
end

M = Q*M;
beta_u = betahat + M*x'*(y - exp(prob)./(exp(prob)+1))/n;
beta_u = c * beta_u;

V = diag(c*M*Sigma*M'*c);
%V = diag(M*Sigma*M');

width = norminv(1-(1-level)/2)*sqrt(V/n);
CI_u = beta_u + width;
CI_l = beta_u - width;


