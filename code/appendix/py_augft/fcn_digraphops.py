import numpy as np

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
    U, Lmd = np.linalg.eig(D - Ap)
    idxp = np.argsort(np.real(Lmd))
    U = U[:, idxp]
    Lmd = np.diag(np.real(Lmd[idxp]))
    s = np.sign(np.mean(U[:, 0]))
    U = s * U

    # 拡張GFT基底（交代成分）と拡張グラフ周波数（交代成分）
    N = Am.shape[0]
    isOdd = N % 2
    V, Gma = np.linalg.eig(-Am)
    gma = np.imag(Gma.diagonal())
    idxgma = np.argsort(gma)[::-1]
    idxgmap = idxgma[:N // 2]
    idxgmam = idxgma[N // 2 + 1 + isOdd:][::-1]
    idxsrtdgma = np.ravel(np.vstack((idxgmap, idxgmam)))
    Gma = 1j * np.diag(gma[idxsrtdgma])
    V = V[:, idxsrtdgma]
    r = np.linalg.matrix_rank(Am)
    Vnz = V[:, :r]
    Qnz = np.reshape(np.array([[1, 1], [1j, -1j]]) @ np.reshape(Vnz.T, (2, -1)) / np.sqrt(2), (N, -1)).T
    Q = Qnz
    Gma_ = Gma[:r, :r]
    Sgm = np.reshape(np.array([[0, 1j], [1j, 0]]) @ np.reshape(Gma_.T, (2, -1)), (r, -1)).T

    # 確認
    assert np.linalg.norm(U @ U.T - np.eye(U.shape[0])) / U.size < 1e-3, "Uの直交性"
    assert np.linalg.norm(Q.T @ Q - np.eye(Q.shape[1])) / Q.size < 1e-3, "Qの直交性"
    assert np.linalg.norm((D - Ap) - U @ Lmd @ U.T) / D.size < 1e-3, "(D-Ap)の分解"
    assert np.linalg.norm(Am + Q @ Sgm @ Q.T) / Am.size < 1e-3, "-Amの分解"

    return U, Q, C, D, L, Lmd, Sgm

"""
function [U,Q,C,D,L,Lmd,Sgm] = fcn_digraphops(A)
% 隣接行列からの拡張グラフ作用素生成
%
% 入力
%
%   A: 隣接行列
%
% 出力
% 
%   U: 拡張GFT 基底（対称成分）
%   Q: 拡張GFT 基底（交代成分）
%   C: 準隣接行列
%   D: 度数行列
%   L: 準グラフラプラシアン
%   Lmd: 拡張グラフ周波数（対称成分）
%   Sgm: 拡張グラフ周波数（交代成分）
%

arguments
    A (:,:) double
end

%　引数がグラフの場合，隣接行列を抽出
if isa(A,'digraph') || isa(A,'graph')
    A = adjacency(A);
end
if issparse(A)
    A = full(A);
end

% 対称成分
Ap = (A + A.')/2;
% 交代成分
Am = (A - A.')/2;
% 準隣接行列
C = Ap + 1j*Am;
% 度数行列
D = diag(sum(abs(C),2));
% 準グラフララプラシアン
L = D - C;

% 拡張GFT基底（対称成分）と拡張グラフ周波数（対称成分）
[U,Lmd] = eig(D-Ap);
[lmd,idxp] = sort(diag(Lmd),'ascend');
U = U(:,idxp);
Lmd = diag(lmd);
s = sign(mean(U(:,1),'all'));
U = s*U;

% 拡張GFT基底（交代成分）と拡張グラフ周波数（交代成分）
N = size(Am,1);
isOdd = mod(N,2);
[V,Gma] = eig(-Am); % マイナスに修正 2023/2/12
gma = diag(imag(Gma)); % 固有値の虚部を抽出）
[~,idxgma] = sort(gma,'descend'); % 大きい順にソート
idxgmap = idxgma(1:floor(end/2)); % 正の固有値
idxgmam = idxgma(end:-1:floor(end/2)+1+isOdd); % 負の固有値
idxsrtdgma = reshape(vertcat(idxgmap.',idxgmam.'),[],1); % 正負固有値の並べ替え
Gma = 1j*diag(gma(idxsrtdgma)); % 虚数に回復
V = V(:,idxsrtdgma); % 有意な固有ベクトル
r = rank(Am); % 交代行列のランク
Vnz = V(:,1:r); % 非ゼロ固有値に対応する固有ベクトル
Qnz = reshape([1 1; 1j -1j]*reshape(Vnz',2,[])/sqrt(2),[],N).'; 
Q = Qnz; % 変換行列
Gma_ = Gma(1:r,1:r);
Sgm = reshape([0 1j; 1j 0]*reshape(Gma_',2,[]),[],r).'; % ブロック対角化

% 確認
assert(norm(U*U.'-eye(size(U)),'fro')/numel(U)<1e-3,'Uの直交性') % Uの直交性
assert(isempty(Q) || norm(Q.'*Q-eye(size(Q,2)),'fro')/numel(Q)<1e-3,'Qの直交性') % Qの直交性
assert(norm((D-Ap)-U*Lmd*U.','fro')/numel(D)<1e-3,'(D-Ap)の分解') % (D-Ap)の分解
assert(norm(Am+Q*Sgm*Q.','fro')/numel(Am)<1e-3,'-Amの分解') % -Am の分解 2023/2/12 修正
end
"""