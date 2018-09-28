clear
clc
close all

% 主程序
% 20180309 v3 重新改回1小时，可以再进一步加大预测误差，风光的标准差分开，数据更贴近真实

global period off_grid EH1 EH2 EH3 Grid1 IESNUMBER
period = 60 / 60; % 分母是时间间隔

if period == 1
    load '../tmp1.mat' % 有177家公司
    load '../tmp2.mat'
    load '../renewableName.mat'
    load '../solarValue.mat'
    load '../windValue.mat'
elseif period == 4 % 只有174家公司了
    load data_loadValue_15min.mat
    load data_loadName_15min.mat
    load renewableName_15min.mat
    load solarValue_15min.mat
    load windValue_15min.mat
end
load '../gridPriceRecord'
Le_max = [2 , 2.3 , 1.5] * 1000;
Le_dr_rate = [0.1, 0.2, 0.4];
Lh_max = [1 , 2 , 2] * 1000;
Lh_dr_rate = [0.2, 0.3, 0.3];
solar_max = [1 , 0.1 , 0.5] * 1000;
wind_max = [0.3 , 1 , 0.1] * 1000;
% IES1 工业区，电热负荷都比较平，白天稍高，热大于电，负荷型
EH1_Le = loadValue(:,94); % 对应134号公司 %/10
EH1_Lh = loadValue(:,143); % 对应223号公司 %/5
EH1_solarP = solarValue(:,3) ; %*1
EH1_windP = windValue(:,3); %*2

% IES2 商务区，电热负荷白天高，晚上很低，热电相当，电有驼峰，负荷型
EH2_Le = loadValue(:,88) ; % 对应127号公司 %/8
EH2_Lh = loadValue(:,91) ; % 对应130号公司 %/6
EH2_solarP = solarValue(:,1) ;
%         EH2_windP = windValue(:,1);
EH2_windP = [windValue(12:end,1);windValue(1:11,1)];
% IES3 住宅区，电热负荷白天低，晚上高，热电相当，资源丰富型
if period == 1
    EH3_Le = loadValue(:,172) ; % 对应273号公司 %/5
    EH3_Lh = loadValue(:,174) ; % 对应275号公司 %/6
elseif period == 4
    EH3_Le = loadValue(:,169) ; % 对应273号公司
    EH3_Lh = loadValue(:,171) ; % 对应275号公司
end
EH3_solarP = solarValue(:,27) ;
EH3_windP = windValue(:,27);

IESNUMBER = 3; % IES个数

for IES_no = 1 : IESNUMBER
    eval(['EH',num2str(IES_no),'_Le = EH',num2str(IES_no),'_Le / max(EH',num2str(IES_no),'_Le) * Le_max(IES_no) * (1 -Le_dr_rate(IES_no));']);
    eval(['EH',num2str(IES_no),'_Lh = EH',num2str(IES_no),'_Lh / max(EH',num2str(IES_no),'_Lh) * Lh_max(IES_no) * (1 -Lh_dr_rate(IES_no));']);
    eval(['EH',num2str(IES_no),'_solarP = EH',num2str(IES_no),'_solarP / max(EH',num2str(IES_no),'_solarP) * solar_max(IES_no);']);
    eval(['EH',num2str(IES_no),'_windP = EH',num2str(IES_no),'_windP / max(EH',num2str(IES_no),'_windP) * wind_max(IES_no);']);
    eval(['EH',num2str(IES_no),'_solarP_rate = solar_max(IES_no);']);
    eval(['EH',num2str(IES_no),'_windP_rate = wind_max(IES_no);']);
end

% 简要画图，判断净负荷，包括电与热

EH1_Le_jing = EH1_Le-EH1_solarP-EH1_windP;
EH2_Le_jing = EH2_Le-EH2_solarP-EH2_windP;
EH3_Le_jing = EH3_Le-EH3_solarP-EH3_windP;

%{
        for IES_no = 1 : IESNUMBER
            subplot(3,2,IES_no * 2 - 1);hold on
            eval(['plot(EH', num2str(IES_no) , '_Le)']);
            eval(['plot(EH',num2str(IES_no) , '_Lh)']);
            eval(['plot(EH',num2str(IES_no) , '_Le_jing)']);
            legend('电负荷','热负荷','净电负荷');
            subplot(3,2, IES_no * 2);hold on
            eval(['plot(EH', num2str(IES_no), '_Le_jing/0.35-EH', num2str(IES_no),'_Lh/0.45);']);
        end
%}
clear loadName loadValue renewableName solarValue windValue


global minMarketPrice maxMarketPrice priceNumbers step

minMarketPrice = 0.1;
maxMarketPrice = 1;
step = 0.1; %只有在单次出清的时候用得到
priceNumbers = (maxMarketPrice - minMarketPrice)/step + 1; %一个投标向量的长度（点数），段数 + 1 = 点数
marketInfo = [minMarketPrice; maxMarketPrice; step; priceNumbers];

% 电网特性
global elePrice
%分时电价
%{
        elePrice = ones(24*period,1) .* 0.6268;
        elePrice(0+1 : 8*period) = ones(8*period,1) .* 0.3089;
        elePrice(8*period+1 : 12*period) = ones(4*period,1) .* 1.0447;
        elePrice(17*period+1 : 21*period) = ones(4*period,1) .* 1.0447;
%}
%实时电价
elePrice = gridPriceRecord( 49 : 49 + 24 - 1)';
elePrice = ( elePrice - min(elePrice) ) / (max(elePrice) - min(elePrice)) * 0.8 + 0.2 ;
clear gridPriceRecord

% 气网特性
global gasPrice1 gasPrice3  % 不需要全局变量吧
gasPrice1 = 0.334; % 3.3元每立方米换算后的值
gasPrice3 = 0.284; % 2.8元每立方米换算后的值
gasLimit1 = 1e6; %暂时不考虑回售天然气
gasLimit2 = 1e6;
gasLimit3 = 1e6;
% gasLimit_total = 150; %暂时无法对天然气总量做单独约束，只能默认是各支线约束的和

%CHP的参数
CHP1_para = [0.30, 0.42, 1400, 0]; % CHP_GE_eff_in, CHP_GH_eff_in, CHP_Prate_in, CHP_Pmin_Rate_in
CHP2_para = [0.35, 0.45, 1, 0];
CHP3_para = [0.28, 0.56, 1200, 0];

%锅炉
Boiler1_para = [0.90; Lh_max(1)/0.9]; % Boiler_eff_in, Boiler_Prate_in
Boiler2_para = [0.90; Lh_max(2)/0.9];
Boiler3_para = [0.90; Lh_max(3)/0.9];

%电储能和热储能
% ES_totalC_in, ES_maxSOC_in, ES_minSOC_in, ES_currentSOC_in, ES_targetSOC_in, ES_chargeTime, ES_eff_in
%{
        ES1_para = [2000, 0.8, 0.2, 0.4, 0.4, 6, 0.9];
        ES2_para = [500, 0.85, 0.15, 0.4, 0.4, 6, 0.9];
        ES3_para = [10, 0.85, 0.15, 0.4, 0.4, 6, 0.9];
        % HS_totalC_in, HS_maxSOC_in, HS_minSOC_in, HS_currentSOC_in, HS_targetSOC_in, HS_chargeTime, HS_eff_in
        HS1_para = [2000, 0.9, 0.1, 0.5, 0.5, 5, 0.9];
        HS2_para = [1500, 0.9, 0.1, 0.5, 0.5, 5, 0.9];
        HS3_para = [10, 0.9, 0.1, 0.5, 0.5, 5, 0.9];
%}
ES1_para = [1, 0.85, 0.15, 0.2, 0.4, 6, 0.9];
ES2_para = [1000, 0.85, 0.15, 0.2, 0.4, 6, 0.9];
ES3_para = [2000, 0.85, 0.15, 0.2, 0.4, 6, 0.9];
% HS_totalC_in, HS_maxSOC_in, HS_minSOC_in, HS_currentSOC_in, HS_targetSOC_in, HS_chargeTime, HS_eff_in
HS1_para = [2000, 0.9, 0.1, 0.6, 0.5, 5, 0.9];
HS2_para = [1000, 0.9, 0.1, 0.6, 0.5, 5, 0.9];
HS3_para = [1, 0.9, 0.1, 0.5, 0.6, 0.5, 0.9];

% 负荷和风光预测误差
dev_L = 3/100; %百分数 1
dev_PV = 10/100; %5
dev_WT = 15/100;
% seedNumber = 0;
% rand 生成均匀分布的伪随机数 分布在（0~1）之间
% randn 生成标准正态分布的伪随机数 （均值为0，标准差为1）（乘以的系数是标准差，不是方差！）
randn('seed', 10);

%电
%各IES各时段最大可平移负荷功率
EH1_Le_drP_rate = Le_max(1) * Le_dr_rate(1);
EH2_Le_drP_rate = Le_max(2) * Le_dr_rate(2);
EH3_Le_drP_rate = Le_max(3) * Le_dr_rate(3);
%各IES总可平移负荷电量
EH1_Le_drP_total = EH1_Le_drP_rate * 10;
EH2_Le_drP_total = EH2_Le_drP_rate * 10;
EH3_Le_drP_total = EH3_Le_drP_rate * 10;

%热
%各IES各时段最大可平移负荷功率
EH1_Lh_drP_rate = Lh_max(1) * Lh_dr_rate(1);
EH2_Lh_drP_rate = Lh_max(2) * Lh_dr_rate(2);
EH3_Lh_drP_rate = Lh_max(3) * Lh_dr_rate(3);
%各IES总可平移负荷电量
EH1_Lh_drP_total = EH1_Lh_drP_rate * 10;
EH2_Lh_drP_total = EH2_Lh_drP_rate * 10;
EH3_Lh_drP_total = EH3_Lh_drP_rate * 10;

singleLimit = Le_max * 1.2;
totalLimit = mean(EH1_Le_jing) + mean(EH2_Le_jing) + mean(EH3_Le_jing)+...
    (EH1_Le_drP_total +EH1_Le_drP_total + EH1_Le_drP_total )/(24 * period) ;
reverseRate = 4;
% 支线: 下级向上级购电、售电约束，再加一个线损率5-7%
eleLimit1 = [singleLimit(1), -singleLimit(1)/reverseRate, 0.94];
eleLimit2 = [singleLimit(2), -singleLimit(2)/reverseRate, 0.94];
eleLimit3 = [singleLimit(3), -singleLimit(3)/reverseRate, 0.94];
eleLimit_total = [totalLimit * 1.1, -totalLimit/reverseRate]; % 馈线

% 用于存放最终优化结果
result_ES_SOC = zeros(24 * period + 1 , IESNUMBER);
result_HS_SOC = zeros(24 * period + 1 , IESNUMBER);
result_Ele = zeros(24 * period , IESNUMBER);
result_CHP_G = zeros(24 * period , IESNUMBER);
result_Boiler_G = zeros(24 * period , IESNUMBER);
result_ES_discharge = zeros(24 * period , IESNUMBER);
result_ES_charge = zeros(24 * period , IESNUMBER);
result_HS_discharge = zeros(24 * period , IESNUMBER);
result_HS_charge = zeros(24 * period , IESNUMBER);

% 电网
Grid1 = Grid_171118(eleLimit_total);
% EH参数与实例化
EH1 = EH_local_170828_v3(eleLimit1, gasLimit1, EH1_Le, EH1_Lh, EH1_solarP, EH1_windP, CHP1_para, Boiler1_para, ES1_para, HS1_para, dev_L, dev_PV, dev_WT, EH1_solarP_rate, EH1_windP_rate, EH1_Le_drP_rate, EH1_Le_drP_total, EH1_Lh_drP_rate, EH1_Lh_drP_total);
EH2 = EH_local_170828_v3(eleLimit2, gasLimit2, EH2_Le, EH2_Lh, EH2_solarP, EH2_windP, CHP2_para, Boiler2_para, ES2_para, HS2_para, dev_L, dev_PV, dev_WT, EH2_solarP_rate, EH2_windP_rate, EH2_Le_drP_rate, EH2_Le_drP_total, EH2_Lh_drP_rate, EH2_Lh_drP_total);
EH3 = EH_local_170828_v3(eleLimit3, gasLimit3, EH3_Le, EH3_Lh, EH3_solarP, EH3_windP, CHP3_para, Boiler3_para, ES3_para, HS3_para, dev_L, dev_PV, dev_WT, EH3_solarP_rate, EH3_windP_rate, EH3_Le_drP_rate, EH3_Le_drP_total, EH3_Lh_drP_rate, EH3_Lh_drP_total);

