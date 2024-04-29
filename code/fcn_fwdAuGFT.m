function X = fcn_fwdAuGFT(x,U,Q,alpha) %beta)

if nargin < 4
    %beta = 1/sqrt(2);
    alpha = 1-1/sqrt(2);
end

% Forward AuGFT
I = eye(size(U));
%alpha = (1-sqrt(1-beta^2));
beta = sqrt(alpha*(2-alpha));
fwdAuGFT = [(I-alpha*(Q*Q.'))*U beta*Q].';
X = fwdAuGFT*x;
end