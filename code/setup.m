% Setup script
%
% Requirements: MATLAB R2022b
%
% Copyright (c) 2022, Hotaka KITAMURA and Shogo MURAMATSU
%
% All rights reserved.
%
% Contact address: Shogo MURAMATSU,
%    Faculty of Engineering, Niigata University,
%    8050 2-no-cho Ikarashi, Nishi-ku,
%    Niigata, 950-2181, JAPAN
%
% http://msiplab.eng.niigata-u.ac.jp/

%% Create RESULTS folder
if ~exist("../results","dir")
    mkdir("../results")
else
    disp("../results exists.")
end

%% Download GSPBox package
GSPBOX_VERSION = "0.7.5";
GSPBOX_TAG = "0.75";
GSPBOX_ROOT = "gspbox";
if ~exist(GSPBOX_ROOT,"dir")
    unzip("https://github.com/epfl-lts2/gspbox/releases/download/"+...
        GSPBOX_TAG +...
        "/gspbox-" +GSPBOX_VERSION+".zip")
end
addpath(GSPBOX_ROOT)
CURRENT_DIR = pwd;
cd(GSPBOX_ROOT)
gsp_start
gsp_install
cd(CURRENT_DIR)

%% Download Global Map Japan 
GMAPJP_VERSION = "2.2";
GMAPJP_DATA = "gm-jpn-all_u_" + GMAPJP_VERSION.replace(".","_");
GMAPJP_ZIP = "https://www1.gsi.go.jp/geowww/globalmap-gsi/download/data/gm-japan/"+GMAPJP_DATA+".zip";
CURRENT_DIR = pwd;
DATA_DIR = "../data";
if ~exist("../data/"+GMAPJP_DATA,"dir")
    cd(DATA_DIR)
    disp("Downloading " + GMAPJP_ZIP)
    unzip(GMAPJP_ZIP,".")
    cd(CURRENT_DIR)
else
    disp("../data/"+GMAPJP_DATA+" exists.")
end

