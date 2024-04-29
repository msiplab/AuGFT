function x = fcn_invAuGFT(X,U,Q,alpha) %beta)

if nargin < 4
    % beta = 1/sqrt(2);
    alpha = 1-1/sqrt(2);
end

% Inverse AuGFT
I = eye(size(U));
%alpha = (1-sqrt(1-beta^2));
beta = sqrt(alpha*(2-alpha));
invAuGFT = [(I-alpha*(Q*Q.'))*U beta*Q];
x = invAuGFT*X;
end