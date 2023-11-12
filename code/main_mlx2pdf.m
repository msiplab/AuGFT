function main_mlx2pdf(varargin)
%% ライブスクリプトからPDFへの変換
%
% Copyright (c) Shogo MURAMATSU, 2023
% All rights resereved

srcDir = fullfile(pwd,'./');
dstDir = './';
isVerbose = true;

%% ファイルの取得
if nargin < 1
    list = ls([srcDir '*.mlx']);
else
    list = ls(sprintf('%s/%s*.mlx',srcDir,varargin{1}));    
end

%% ファイルの変換
for idx = 1:size(list,1)
    % ファイル名の抽出
    [~,fname,~] = fileparts(list(idx,:));
    % スクリプトへ変換
    mlx2pdf(srcDir,fname,dstDir,isVerbose)
    % 変換後のスクリプトの内容
    %open(fullfile(dstDir,[fname '.m']))
end
end

function mlx2pdf(srcDir,fname,dstDir,isVerbose)
%M2MLX
%
% Copyright (c) Shogo MURAMATSU, 2020
% All rights reserved

% 出力フォルダを準備
if exist(dstDir,'dir') ~= 7
    mkdir(dstDir)
end

% ライブスクリプトからスクリプトに変換
if exist(dstDir,'dir') == 7
    srcFile = fullfile(srcDir, [fname '.mlx']);
    dstFile = fullfile(dstDir, [fname '.pdf']);
    if exist(srcFile,'file') == 2
        matlab.internal.liveeditor.openAndConvert(srcFile, dstFile);
        if isVerbose
            fprintf('Saved %s as PDF\n',srcFile);
        end
    else
        me = MException('MsipException:NoSuchFiler', ...
        sprintf('%s does not exist',strrep(srcFile,'\','\\')));
        throw(me)
    end
else
    me = MException('MsipException:NoSuchFolder', ...
    sprintf('%s does not exist',strrep(dstDir,'\','\\')));
    throw(me)
end
end