function [U,F,C,D,L,Lmd] = fcn_paragraphops(A)
% 隣接行列からの準グラフ作用素生成
%
% 入力
%   A: 隣接行列
%
% 出力
%
%   U: ParaGFT 基底 (ParaGFT U'*x)
%   F: 実フィルタ制約行列 ( bar h = F h )
%   C: 準隣接行列
%   D: 度数行列
%   L: 準グラフラプラシアン
%   Lmd: 準グラフ周波数
%

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
% ParaGFT基底と準グラフ周波数
[U,Lmd] = eig(L);
[lmd,idxp] = sort(diag(Lmd),'ascend');
U = U(:,idxp);
Lmd = diag(lmd);
% 実フィルタ制約行列
Xi = zeros(size(U));
for jRow = 1:size(Xi,1)
    u_j = U(:,jRow);
    for iCol = 1:size(Xi,2)
        u_i = U(:,iCol);
        Xi(jRow,iCol) = u_j'*conj(u_i);
    end
end
F = abs(Xi).^2;
%
% 確認
assert(norm(U*U'-eye(size(U)),'fro')/numel(U)<1e-3,'Uのユニタリ性') % Uのユニタリ性
end