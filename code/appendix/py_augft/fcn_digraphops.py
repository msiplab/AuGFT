import numpy as np
import networkx as nx

def fcn_digraphops(A):
    """
    隣接行列からの拡張グラフ作用素生成

    入力

        A: 隣接行列

    出力

        U: 拡張GFT 基底（対称成分）
        Q: 拡張GFT 基底（交代成分）
        C: 準隣接行列
        D: 度数行列
        L: 準グラフラプラシアン
        Lmd: 拡張グラフ周波数（対称成分）
        Sgm: 拡張グラフ周波数（交代成分）
    """

    # 引数がグラフの場合，隣接行列を抽出
    if isinstance(A, np.ndarray):
        A = A.copy()
    elif isinstance(A, np.matrix):
        A = np.array(A)
    elif isinstance(A, (nx.Graph, nx.DiGraph)):
        A = nx.adjacency_matrix(A).toarray()
    else:
        raise ValueError("Invalid input type for A")

    # 対称成分
    Ap = (A + A.T) / 2
    # 交代成分
    Am = (A - A.T) / 2
    # 準隣接行列
    C = Ap + 1j * Am
    # 度数行列
    D = np.diag(np.sum(np.abs(C), axis=1))
    # 準グラフラプラシアン
    L = D - C

    # 拡張GFT基底（対称成分）と拡張グラフ周波数（対称成分）
    Lmd, U = np.linalg.eigh(D - Ap)
    idxp = np.argsort(np.real(Lmd))
    U = U[:, idxp]
    Lmd = np.diag(np.real(Lmd[idxp]))
    s = np.sign(np.mean(U[:, 0]))
    U = s * U

    # 拡張GFT基底（交代成分）と拡張グラフ周波数（交代成分）
    N = Am.shape[0]
    isOdd = N % 2
    #
    Gma, V = np.linalg.eigh(1j*Am) # np.linalg.eig(-Am)
    gma = np.real(Gma) # np.imag(Gma) # 固有値の虚部を抽出    
    idxgma = np.argsort(gma)[::-1] # 大きい順にソート
    idxgmap = idxgma[:N // 2] # 正の固有値
    idxgmam = idxgma[N // 2 + isOdd:][::-1] # 負の固有値
    idxsrtdgma = np.empty((idxgmap.size + idxgmam.size,), dtype=idxgmap.dtype)
    idxsrtdgma[0::2] = idxgmap
    idxsrtdgma[1::2] = idxgmam

    Gma = 1j * np.diag(gma[idxsrtdgma]) # 虚部に回復
    V = V[:, idxsrtdgma] # 有意な固有ベクトル
    r = np.linalg.matrix_rank(Am) # 交代行列のランク
    Vnz = V[:, :r] # 非ゼロ固有値に対する固有ベクトル
    Qnz = np.real(Vnz.reshape((-1,2)) @ np.array([ [-1, 1j], [-1, -1j] ])).reshape((N,-1)) / np.sqrt(2)
    Q = Qnz # 変換行列
    Gma_ = Gma[:r, :r] 
    Sgm = np.real(Gma_.reshape((-1,2)) @ np.array([ [0, -1j], [-1j, 0] ])).reshape(r,-1) # ブロック対角化

    # 確認
    assert np.linalg.norm(U @ U.T - np.eye(U.shape[0])) / U.size < 1e-3, "Uの直交性"
    assert np.linalg.norm(Q.T @ Q - np.eye(Q.shape[1])) / Q.size < 1e-3, "Qの直交性"
    assert np.linalg.norm((D - Ap) - U @ Lmd @ U.T) / D.size < 1e-3, "(D-Ap)の分解"
    assert np.linalg.norm(Am + Q @ Sgm @ Q.T) / Am.size < 1e-3, "-Amの分解"

    return U, Q, C, D, L, Lmd, Sgm
