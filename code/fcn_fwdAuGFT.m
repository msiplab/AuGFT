function [X,fwdAuGFT] = fcn_fwdAuGFT(x,U,Q,alpha) 

if nargin < 4
    alpha = 1-1/sqrt(2);
end

% Forward AuGFT
I = eye(size(U));
beta = sqrt(alpha*(2-alpha));
fwdAuGFT = [(I-alpha*(Q*Q.'))*U beta*Q].';
X = fwdAuGFT*x;
end