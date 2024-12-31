import numpy as np

def fcn_invAuGFT(X, U, Q, alpha=(1-1/np.sqrt(2))):
    
    # Inverse AuGFT
    I = np.eye(U.shape[0])
    beta = np.sqrt(alpha * (2 - alpha))
    invAuGFT = np.hstack(((I - alpha * (Q@Q.T)) @ U, beta*Q))
    x = invAuGFT@X

    return (x, invAuGFT)