function X = fcn_fwdAuGFT(x,U,Q,beta)

if nargin < 4
    beta = 1/sqrt(2);
end

% Forward AuGFT
I = eye(size(U));
alpha = (1-sqrt(1-beta^2));
fwdAuGFT = [(I-alpha*(Q*Q.'))*U beta*Q].';
X = fwdAuGFT*x;
end