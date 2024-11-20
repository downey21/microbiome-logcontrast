function [betahat]= bycvx(penalize, y,x, constr2, lambda)
tic

[n p] = size(x);

if penalize == 0
    if (sum(abs(constr2))==0)
        cvx_begin quiet
        cvx_precision low
        variable betahat(p)
            minimize( -(y'*x*betahat - sum(log(1+exp(x*betahat))))+n*lambda*norm(betahat,1) )
        cvx_end
    else
        cvx_begin quiet
        cvx_precision low
        variable betahat(p)
            minimize( -(y'*x*betahat - sum(log(1+exp(x*betahat))))+n*lambda*norm(betahat,1) )
            subject to
            constr2*betahat == 0
        cvx_end
    end
else
    if (sum(abs(constr2))==0)
        cvx_begin quiet
        cvx_precision low
        variable betahat(p)
            minimize( -(y'*x*betahat - sum(log(1+exp(x*betahat))))+n*lambda*norm(betahat(2:p),1) )
        cvx_end
    else
        cvx_begin quiet
        cvx_precision low
        variable betahat(p)
            minimize( -(y'*x*betahat - sum(log(1+exp(x*betahat))))+n*lambda*norm(betahat(2:p),1) )
            subject to
            constr2*betahat == 0
        cvx_end
    end
end

toc















