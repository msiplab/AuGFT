function [x,invAuGFT] = fcn_invAuGFT(X,U,Q,alpha) 

if nargin < 4
    alpha = 1-1/sqrt(2);
end

% Inverse AuGFT
I = eye(size(U));
beta = sqrt(alpha*(2-alpha));
invAuGFT = [(I-alpha*(Q*Q.'))*U beta*Q];
x = invAuGFT*X;
end