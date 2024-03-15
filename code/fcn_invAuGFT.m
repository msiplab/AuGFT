function x = fcn_invAuGFT(X,U,Q,beta)

if nargin < 4
    beta = 1/sqrt(2);
end

% Inverse AuGFT
I = eye(size(U));
alpha = (1-sqrt(1-beta^2));
invAuGFT = [(I-alpha*(Q*Q.'))*U beta*Q];
x = invAuGFT*X;
end