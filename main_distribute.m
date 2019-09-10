clc;clear;
global caseType couldExport
caseType = 2;
couldExport = 1;
para_init;
off_grid = 0; % 0表示正常运行，1表示IES1离网
priceArray = elePrice;
priceArray_record = zeros( 24 * period , 3);
% 2-stage
%日前优化
isDA = 1;
% 次梯度法
all_temporal;

%日内优化
isDA = 0;
%单时段滚动求解
% 二分法
single_temporal;
isCentral = 0;
main_handle_171013_v2

