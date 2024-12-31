import numpy as np

def fcn_invAuGFT(X, U, Q, **kwargs):

    if 'alpha' in kwargs:
        alpha = kwargs.get('alpha')
    else:
        alpha = 1 - 1 / np.sqrt(2)
    
    # Inverse AuGFT
    I = np.eye(U.shape[0])
    beta = np.sqrt(alpha * (2 - alpha))
    invAuGFT = np.concatenate(((I - alpha * np.dot(Q, Q.T)) * U, beta * Q), axis=1)
    x = np.dot(invAuGFT, X)
    return x, invAuGFT