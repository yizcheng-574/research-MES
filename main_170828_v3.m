clear
clc
close all

% 主程序
% 20180309 v3 重新改回1小时，可以再进一步加大预测误差，风光的标准差分开，数据更贴近真实

global period
period = 60 / 15; % 分母是时间间隔

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

% IES1 工业区，电热负荷都比较平，白天稍高，热大于电，负荷型
EH1_Le = loadValue(:,94) .* 4; % 对应134号公司 %/10
EH1_Lh = loadValue(:,143) .* 8; % 对应223号公司 %/5
EH1_solarP = solarValue(:,3) ./1000 .* 50; %*1
EH1_windP = windValue(:,3)./1000 .* 30; %*2
EH1_solarP_rate = 5000;
EH1_windP_rate = 1500;
EH1_solarP(EH1_solarP>EH1_solarP_rate) = EH1_solarP_rate;
EH1_windP(EH1_windP>EH1_windP_rate) = EH1_windP_rate;

% IES2 商务区，电热负荷白天高，晚上很低，热电相当，电有驼峰，负荷型
EH2_Le = loadValue(:,88) ./ 2.25; % 对应127号公司 %/8
EH2_Lh = loadValue(:,91) ./ 2.25*1.4 ; % 对应130号公司 %/6
EH2_solarP = solarValue(:,1) ./1000 .* 2; %*1
EH2_windP = windValue(:,1)./1000 .* 15; %*3
EH2_solarP_rate = 250;
EH2_windP_rate = 500;
EH2_solarP(EH2_solarP>EH2_solarP_rate) = EH2_solarP_rate;
EH2_windP(EH2_windP>EH2_windP_rate) = EH2_windP_rate;

% IES3 住宅区，电热负荷白天低，晚上高，热电相当，资源丰富型
if period == 1
    EH3_Le = loadValue(:,172) .* 3; % 对应273号公司 %/5
    EH3_Lh = loadValue(:,174) .* 2.5; % 对应275号公司 %/6
elseif period == 4
    EH3_Le = loadValue(:,169) ./ 5; % 对应273号公司
    EH3_Lh = loadValue(:,171) ./ 6; % 对应275号公司
end
EH3_solarP = solarValue(:,27) ./1000 .* 9; %*2
EH3_windP = windValue(:,27)./1000 .* 6; %*6
EH3_solarP_rate = 1000;
EH3_windP_rate = 300;
EH3_solarP(EH3_solarP>EH3_solarP_rate) = EH3_solarP_rate;
EH3_windP(EH3_windP>EH3_windP_rate) = EH3_windP_rate;


%{
% 简要画图，判断净负荷，包括电与热
EH1_Le_jing = EH1_Le-EH1_solarP-EH1_windP;
EH2_Le_jing = EH2_Le-EH2_solarP-EH2_windP;
EH3_Le_jing = EH3_Le-EH3_solarP-EH3_windP;
figure
hold on
plot(EH1_Le)
plot(EH1_Lh,'r')
plot(EH1_Le_jing,'k')
figure
hold on
plot(EH2_Le)
plot(EH2_Lh,'r')
plot(EH2_Le_jing,'k')
figure
hold on
plot(EH3_Le)
plot(EH3_Lh,'r')
plot(EH3_Le_jing,'k')
%}
clear loadName loadValue renewableName solarValue windValue


global minMarketPrice maxMarketPrice priceNumbers step

minMarketPrice = 0;
maxMarketPrice = 1.5;
step = 0.1; %只有在单次出清的时候用得到
priceNumbers = (maxMarketPrice - minMarketPrice)/step + 1; %一个投标向量的长度（点数），段数 + 1 = 点数
marketInfo = [minMarketPrice; maxMarketPrice; step; priceNumbers];

% 电网特性
global elePrice
elePrice = ones(24*period,1) .* 0.6268;
elePrice(0+1 : 8*period) = ones(8*period,1) .* 0.3089;
elePrice(8*period+1 : 12*period) = ones(4*period,1) .* 1.0447;
elePrice(17*period+1 : 21*period) = ones(4*period,1) .* 1.0447;

singleLimit = [16000, 800, 2500]; %[16000, 800, 2500] [0 0 0] [16000, 800, 2500]./19300.*16000
totalLimit = 16000; 
reverseRate = 4;
% 支线: 下级向上级购电、售电约束，再加一个线损率5-7%
eleLimit1 = [singleLimit(1), -singleLimit(1)/reverseRate, 0.94];
eleLimit2 = [singleLimit(2), -singleLimit(2)/reverseRate, 0.94];
eleLimit3 = [singleLimit(3), -singleLimit(3)/reverseRate, 0.94];
eleLimit_total = [totalLimit, -totalLimit/reverseRate]; % 馈线

% 气网特性
% global gasPrice1 gasPrice3  % 不需要全局变量吧
gasPrice1 = 0.334; % 3.3元每立方米换算后的值
gasPrice3 = 0.284; % 2.8元每立方米换算后的值
gasLimit1 = 1e6; %暂时不考虑回售天然气
gasLimit2 = 1e6;
gasLimit3 = 1e6;
% gasLimit_total = 150; %暂时无法对天然气总量做单独约束，只能默认是各支线约束的和

%CHP的参数
CHP1_para = [0.35, 0.45, 30000, 0.4]; % CHP_GE_eff_in, CHP_GH_eff_in, CHP_Prate_in, CHP_Pmin_Rate_in
CHP2_para = [0.35, 0.45, 1500, 0.45];
CHP3_para = [0.35, 0.45, 5000, 0.45];

%锅炉
Boiler1_para = [0.90; 1e5]; % Boiler_eff_in, Boiler_Prate_in
Boiler2_para = [0.90; 1e5];
Boiler3_para = [0.90; 1e5];

%电储能和热储能
% ES_totalC_in, ES_maxSOC_in, ES_minSOC_in, ES_currentSOC_in, ES_targetSOC_in, ES_chargeTime, ES_eff_in
% 0.096*200, 0.85, 0.15, 0.5, 0.5, 0.024*200, 0.9
ES1_para = [50000, 0.8, 0.2,      0.4, 0.4,   10, 0.9]; 
ES2_para = [3000, 0.85, 0.15,      0.4, 0.4,   10, 0.9]; 
ES3_para = [10000, 0.85, 0.15,    0.4, 0.4,   10, 0.9]; 
% HS_totalC_in, HS_maxSOC_in, HS_minSOC_in, HS_currentSOC_in, HS_targetSOC_in, HS_chargeTime, HS_eff_in
% 5*4, 0.85, 0.15, 0.5, 0.5, 1.5*4, 0.9
HS1_para = [6000, 0.9, 0.1, 0.5, 0.5, 5, 0.9];
HS2_para = [1000, 0.9, 0.1, 0.5, 0.5, 5, 0.9];
HS3_para = [6000, 0.9, 0.1, 0.5, 0.5, 5, 0.9]; 

% 负荷和风光预测误差
dev_L = 3/100; %百分数 1
dev_PV = 10/100; %5
dev_WT = 15/100;
% seedNumber = 0;
% rand 生成均匀分布的伪随机数 分布在（0~1）之间
% randn 生成标准正态分布的伪随机数 （均值为0，标准差为1）（乘以的系数是标准差，不是方差！）
randn('seed', 10);

% 柴油发电机
% DG2 = DieselGenerator_171123(100, 0.001, 0.125, 0.3);


% 集中式优化：单次
%{
EH1_Le_final = [210.301397758191;230.823676166166;214.108947989989;209.606808045361;196.528718871749;200.527716727458;171.826840261211;215.434468674951;351.421093302298;385.510286550270;408.046130985166;351.424185163977;253.314128307174;312.468981661435;382.019387090551;393.974804830332;395.631704530858;295.363396526140;294.358987619370;322.389351836491;303.152990891480;218.531577322458;250.288679710111;228.503601958543;];
EH1_Lh_final = [575.110768433430;566.343928391538;565.561753170696;558.747459140154;559.049250924452;577.664251332143;573.490976377669;576.759582520053;594.808837746766;572.100185998495;584.937565437560;613.915598291639;591.218527787315;605.178488236794;616.083673626242;630.747911468568;643.183145122572;653.889166802478;643.201515608216;606.371220233789;632.736788941651;644.504932902377;600.004163684607;575.254226710805;];

EH2_Le_final = [410.408897758191;465.816176166166;452.198947989989;420.799308045361;403.801218871749;369.110216727458;332.471840261211;328.449468674951;451.846093302298;567.285286550270;578.236130985166;550.179185163977;440.811628307173;505.541481661435;530.489387090552;491.007304830332;494.156704530858;568.975896526140;543.006487619370;559.366851836491;499.457990891480;512.581577322458;494.963679710112;398.756101958543;];
EH2_Lh_final = [576.424768433430;587.953928391538;602.429753170696;595.615959140154;588.405250924452;562.561751332143;567.056976377669;601.085082520053;639.772837746766;645.834185998495;635.211065437560;629.742598291639;618.848027787315;641.226488236794;624.429173626242;643.937911468568;627.137145122572;620.676666802478;612.978515608216;597.720720233789;595.010788941651;620.424432902377;583.820163684607;582.032226710805;];

[EH1_f, EH1_ub, EH1_lb, EH1_Aeq, EH1_beq, EH1_A, EH1_b, EH1_A_eleLimit_total] = OptMatrix(eleLimit1, gasLimit1, EH1_Le_final, EH1_Lh_final, CHP1_para, Boiler1_para, ES1_para, HS1_para, elePrice, gasPrice);
% [x,fval,exitflag,output,lambda] = linprog(EH1_f, EH1_A, EH1_b, EH1_Aeq, EH1_beq, EH1_lb, EH1_ub) % 单一集中式优化测试
[EH2_f, EH2_ub, EH2_lb, EH2_Aeq, EH2_beq, EH2_A, EH2_b, EH2_A_eleLimit_total] = OptMatrix(eleLimit2, gasLimit2, EH2_Le_final, EH2_Lh_final, CHP2_para, Boiler2_para, ES2_para, HS2_para, elePrice, gasPrice);
% [x,fval,exitflag,output,lambda] = linprog(EH2_f, EH2_A, EH2_b, EH2_Aeq, EH2_beq, EH2_lb, EH2_ub) % 单一集中式优化测试

time = 24; %总时间段
number = 2;
var = time * 7;
totalVar = time * 7 * number; %总变量数
%第1,2,3组time是购电量、CHP购气量、锅炉购气量，第4-7组time是储电、储热的放、充功率

f = [EH1_f; EH2_f];
ub = [EH1_ub; EH2_ub];
lb = [EH1_lb; EH2_lb];
Aeq = [EH1_Aeq, zeros(2,var); zeros(2,var), EH2_Aeq];
beq = [EH1_beq; EH2_beq];

% 不加额外约束
% A = [EH1_A, zeros(var,var); zeros(var,var), EH2_A];
% b = [EH1_b; EH2_b];

% 需要额外增加一个购电量的上、下限约束
A1 = [EH1_A, zeros(var,var); zeros(var,var)  EH2_A];
b1 = [EH1_b; EH2_b];
A2 = [EH1_A_eleLimit_total, EH2_A_eleLimit_total];
b2 = ones(time, 1) .* eleLimit_total(1);
b2_sale = ones(time, 1) .* eleLimit_total(2);
A = [A1; A2; -A2];
b = [b1; b2; -b2_sale];

[x,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub)

x(1:24) + x(24*7+1:24*8)
%}


IESnumber = 3; % IES个数

% 用于存放最终优化结果
result_ES_SOC = zeros(24*period+1,IESnumber);
result_HS_SOC = zeros(24*period+1,IESnumber);
result_Ele = zeros(24*period,IESnumber);
result_CHP_G = zeros(24*period,IESnumber);
result_Boiler_G = zeros(24*period,IESnumber);
result_ES_discharge = zeros(24*period,IESnumber);
result_ES_charge = zeros(24*period,IESnumber);
result_HS_discharge = zeros(24*period,IESnumber);
result_HS_charge = zeros(24*period,IESnumber);



% 集中式优化：滚动 24*period 次
%{
% 初始SOC
result_ES_SOC(1,1) = ES1_para(4);
result_ES_SOC(1,2) = ES2_para(4);
result_ES_SOC(1,3) = ES3_para(4);

result_HS_SOC(1,1) = HS1_para(4);
result_HS_SOC(1,2) = HS2_para(4);
result_HS_SOC(1,3) = HS3_para(4);

for t_current = 1:24*period
    % 更新预测
    [EH1_Le, EH1_Lh, EH1_solarP, EH1_windP] = predict(EH1_Le, EH1_Lh, EH1_solarP, EH1_windP, t_current, dev_L, dev_PV, dev_WT, EH1_solarP_rate, EH1_windP_rate);
    [EH2_Le, EH2_Lh, EH2_solarP, EH2_windP] = predict(EH2_Le, EH2_Lh, EH2_solarP, EH2_windP, t_current, dev_L, dev_PV, dev_WT, EH2_solarP_rate, EH2_windP_rate);
    [EH3_Le, EH3_Lh, EH3_solarP, EH3_windP] = predict(EH3_Le, EH3_Lh, EH3_solarP, EH3_windP, t_current, dev_L, dev_PV, dev_WT, EH3_solarP_rate, EH3_windP_rate);
    % 更新电池SOC
    ES1_para(4) = result_ES_SOC(t_current, 1);
    ES2_para(4) = result_ES_SOC(t_current, 2);
    ES3_para(4) = result_ES_SOC(t_current, 3);
    
    HS1_para(4) = result_HS_SOC(t_current, 1);
    HS2_para(4) = result_HS_SOC(t_current, 2);    
    HS3_para(4) = result_HS_SOC(t_current, 3); 
    
    [EH1_f, EH1_ub, EH1_lb, EH1_A, EH1_b, EH1_A_eleLimit_total] = OptMatrix_rolling_20171010(eleLimit1, gasLimit1, EH1_Le, EH1_Lh, CHP1_para, Boiler1_para, ES1_para, HS1_para, gasPrice1, EH1_windP, EH1_solarP, t_current);
    [EH2_f, EH2_ub, EH2_lb, EH2_A, EH2_b, EH2_A_eleLimit_total] = OptMatrix_rolling_20171010(eleLimit2, gasLimit2, EH2_Le, EH2_Lh, CHP2_para, Boiler2_para, ES2_para, HS2_para, gasPrice1, EH2_windP, EH2_solarP, t_current);
    [EH3_f, EH3_ub, EH3_lb, EH3_A, EH3_b, EH3_A_eleLimit_total] = OptMatrix_rolling_20171010(eleLimit3, gasLimit3, EH3_Le, EH3_Lh, CHP3_para, Boiler3_para, ES3_para, HS3_para, gasPrice3, EH3_windP, EH3_solarP, t_current);

    
    time = 24*period - t_current + 1; %总时间段
    var = time * 7;
    %第1,2,3组time是购电量、CHP购气量、锅炉购气量，第4-7组time是储电、储热的放、充功率
    
    f = [EH1_f; EH2_f; EH3_f];
    ub = [EH1_ub; EH2_ub; EH3_ub];
    lb = [EH1_lb; EH2_lb; EH3_lb];

    % 当前版本没有等式约束了，都是不等式约束
    A1 = [EH1_A, zeros(var+2,var), zeros(var+2,var); 
          zeros(var+2,var), EH2_A, zeros(var+2,var);
          zeros(var+2,var), zeros(var+2,var), EH3_A;];
    b1 = [EH1_b; EH2_b; EH3_b];
    
    % 不增加额外约束
%     A = A1;
%     b = b1;
    
    % 需要额外增加一个对于总线的购电量的上、下限约束
    A2 = [EH1_A_eleLimit_total, EH2_A_eleLimit_total, EH3_A_eleLimit_total];
    b2 = ones(time, 1) .* eleLimit_total(1);
    b2_sale = ones(time, 1) .* eleLimit_total(2);
    A = [A1; A2; -A2];
    b = [b1; b2; -b2_sale];
    
    [x,fval,exitflag,output,lambda] = linprog(f,A,b,[],[],lb,ub);
    if exitflag~=1
        error(['集中式优化失败！时间是',num2str(t_current)])
    end
    
    %测试总线的购电量约束是否满足
    % x(1:24) + x(24*7+1:24*8)
    
    %记录结果
    for i=1 : IESnumber
        % 只执行当前周期的结果
        result_Ele(t_current, i) = x(1 + var*(i-1), 1);
        result_CHP_G(t_current, i) = x(time+1 + var*(i-1), 1);
        result_Boiler_G(t_current, i) = x(time*2+1 + var*(i-1), 1);
        result_ES_discharge(t_current, i) = x(time*3+1 + var*(i-1), 1);
        result_ES_charge(t_current, i) = x(time*4+1 + var*(i-1), 1);
        result_HS_discharge(t_current, i) = x(time*5+1 + var*(i-1), 1);
        result_HS_charge(t_current, i) = x(time*6+1 + var*(i-1), 1);
    end    
    %更新储能状态
    result_ES_SOC(t_current+1, 1) = result_ES_SOC(t_current, 1) - result_ES_discharge(t_current, 1) / ES1_para(7) / ES1_para(1) + result_ES_charge(t_current, 1) * ES1_para(7) / ES1_para(1);
    result_ES_SOC(t_current+1, 2) = result_ES_SOC(t_current, 2) - result_ES_discharge(t_current, 2) / ES2_para(7) / ES2_para(1) + result_ES_charge(t_current, 2) * ES2_para(7) / ES2_para(1);
    result_ES_SOC(t_current+1, 3) = result_ES_SOC(t_current, 3) - result_ES_discharge(t_current, 3) / ES3_para(7) / ES3_para(1) + result_ES_charge(t_current, 3) * ES3_para(7) / ES3_para(1);
    
    result_HS_SOC(t_current+1, 1) = result_HS_SOC(t_current, 1) - result_HS_discharge(t_current, 1) / HS1_para(7) / HS1_para(1) + result_HS_charge(t_current, 1) * HS1_para(7) / HS1_para(1);
    result_HS_SOC(t_current+1, 2) = result_HS_SOC(t_current, 2) - result_HS_discharge(t_current, 2) / HS2_para(7) / HS2_para(1) + result_HS_charge(t_current, 2) * HS2_para(7) / HS2_para(1);
    result_HS_SOC(t_current+1, 3) = result_HS_SOC(t_current, 3) - result_HS_discharge(t_current, 3) / HS3_para(7) / HS3_para(1) + result_HS_charge(t_current, 3) * HS3_para(7) / HS3_para(1);    
end
% 电网
gridClearDemand = - sum(result_Ele,2); %1表示按列求和，2表示按行求和
%}



%分布式方法
global Grid1 EH1 EH2 EH3
% 电网
Grid1 = Grid_171118(eleLimit_total);
% EH参数与实例化
EH1 = EH_local_170828_v3(eleLimit1, gasLimit1, EH1_Le, EH1_Lh, EH1_solarP, EH1_windP, CHP1_para, Boiler1_para, ES1_para, HS1_para, dev_L, dev_PV, dev_WT, EH1_solarP_rate, EH1_windP_rate);
EH2 = EH_local_170828_v3(eleLimit2, gasLimit2, EH2_Le, EH2_Lh, EH2_solarP, EH2_windP, CHP2_para, Boiler2_para, ES2_para, HS2_para, dev_L, dev_PV, dev_WT, EH2_solarP_rate, EH2_windP_rate);
EH3 = EH_local_170828_v3(eleLimit3, gasLimit3, EH3_Le, EH3_Lh, EH3_solarP, EH3_windP, CHP3_para, Boiler3_para, ES3_para, HS3_para, dev_L, dev_PV, dev_WT, EH3_solarP_rate, EH3_windP_rate);


%单次出清
%{
priceArray = elePrice; %由历史数据得到，日前预测电价
priceArray_record = zeros(24*period,2); %记录日前和日内的出清价格

% off_grid = 0; % 0表示正常运行，1表示IES1离网
t_realtime = zeros(24*period,3); %记录每次日内优化所需的时间

for t_current = 1:24*period
    disp(['t_current is ',num2str(t_current)]);
    
    if t_current == 1 %则进行日前优化
        tic
        EH1.predict(0);
        EH2.predict(0);
        EH3.predict(0);
        for pt = t_current : 1 : 24*period
            %投标
            gridDemand = Grid1.zGenerate(elePrice(pt)); % 在当前优化周期内不变
            EH1Demand = EH1.curveGenerate(priceArray, gasPrice1, pt);
            EH2Demand = EH2.curveGenerate(priceArray, gasPrice1, pt);
            EH3Demand = EH3.curveGenerate(priceArray, gasPrice3, pt);
            %聚合
            demand_sum = gridDemand + EH1Demand + EH2Demand + EH3Demand;
            %出清
            priceArray(pt) = clearing(demand_sum, 0); %市场出清得到出清价格，并更新预测电价序列
            %响应出清价格
            Grid1.getClearDemand(priceArray(pt), pt);
            EH1.conditionHandlePrice(priceArray, gasPrice1, pt); %EH收到出清价格，本地再优化一次，更新自身状态
            EH2.conditionHandlePrice(priceArray, gasPrice1, pt);
            EH3.conditionHandlePrice(priceArray, gasPrice3, pt);
        end
        
        priceArray_record(:,1) = priceArray; %记录日前出清价格
        t_dayahead = toc; %记录日前优化所需的总时间（这个时间不准确，因为3个IES实际上应并行计算，而非串行）
        
        % 日前优化的结果
        gridClearDemand = Grid1.getResult;
        [result_Ele(:,1), result_CHP_G(:,1), result_Boiler_G(:,1), result_ES_discharge(:,1), result_ES_charge(:,1), result_HS_discharge(:,1), result_HS_charge(:,1), result_ES_SOC(:,1), result_HS_SOC(:,1), EH1_Le, EH1_Lh, EH1_solarP, EH1_windP] = EH1.getResult;
        [result_Ele(:,2), result_CHP_G(:,2), result_Boiler_G(:,2), result_ES_discharge(:,2), result_ES_charge(:,2), result_HS_discharge(:,2), result_HS_charge(:,2), result_ES_SOC(:,2), result_HS_SOC(:,2), EH2_Le, EH2_Lh, EH2_solarP, EH2_windP] = EH2.getResult;
        [result_Ele(:,3), result_CHP_G(:,3), result_Boiler_G(:,3), result_ES_discharge(:,3), result_ES_charge(:,3), result_HS_discharge(:,3), result_HS_charge(:,3), result_ES_SOC(:,3), result_HS_SOC(:,3), EH3_Le, EH3_Lh, EH3_solarP, EH3_windP] = EH3.getResult;   
    end
    
   
    
    % 日内优化（分为正常运行、有离网）  
%     if off_grid == 0 % 正常运行
        EH1.predict(t_current);
        EH2.predict(t_current);
        EH3.predict(t_current);
        %投标
        gridDemand = Grid1.zGenerate(elePrice(t_current)); % 在当前优化周期内不变
        tic
        EH1Demand = EH1.curveGenerate(priceArray, gasPrice1, t_current);
        t_realtime(t_current,1) = toc;
        tic
        EH2Demand = EH2.curveGenerate(priceArray, gasPrice1, t_current);
        t_realtime(t_current,2) = toc;
        tic
        EH3Demand = EH3.curveGenerate(priceArray, gasPrice3, t_current);
        t_realtime(t_current,3) = toc;
        %聚合
        demand_sum = gridDemand + EH1Demand + EH2Demand + EH3Demand;
        %出清
        priceArray(t_current) = clearing(demand_sum, 0); %更新为实际电价
        %响应出清价格
        Grid1.getClearDemand(priceArray(t_current), t_current);
        EH1.conditionHandlePrice(priceArray, gasPrice1, t_current);
        EH2.conditionHandlePrice(priceArray, gasPrice1, t_current);
        EH3.conditionHandlePrice(priceArray, gasPrice3, t_current);
    
%     else % IES1离网    
%     end
end
priceArray_record(:,2) = priceArray; %记录日内出清价格

% 日内优化的结果
gridClearDemand = Grid1.getResult;
[result_Ele(:,1), result_CHP_G(:,1), result_Boiler_G(:,1), result_ES_discharge(:,1), result_ES_charge(:,1), result_HS_discharge(:,1), result_HS_charge(:,1), result_ES_SOC(:,1), result_HS_SOC(:,1), EH1_Le, EH1_Lh, EH1_solarP, EH1_windP] = EH1.getResult;
[result_Ele(:,2), result_CHP_G(:,2), result_Boiler_G(:,2), result_ES_discharge(:,2), result_ES_charge(:,2), result_HS_discharge(:,2), result_HS_charge(:,2), result_ES_SOC(:,2), result_HS_SOC(:,2), EH2_Le, EH2_Lh, EH2_solarP, EH2_windP] = EH2.getResult;
[result_Ele(:,3), result_CHP_G(:,3), result_Boiler_G(:,3), result_ES_discharge(:,3), result_ES_charge(:,3), result_HS_discharge(:,3), result_HS_charge(:,3), result_ES_SOC(:,3), result_HS_SOC(:,3), EH3_Le, EH3_Lh, EH3_solarP, EH3_windP] = EH3.getResult;
% test内容都放在handle函数中

%}




%20171117 迭代出清，二分法

priceArray = elePrice; %由历史数据得到预测电价 %也是用于迭代的价格变量
priceArray_record = zeros(24*period,2); %一列日前，一列实时
ee = 0.01;

demand_sum = zeros(priceNumbers, 1);
gridClearDemand = zeros(24*period,1);
global iterationNumber
iterationTimes = zeros(24*period, 2); %记录迭代次数

global off_grid
off_grid = 0; % 0表示正常运行，1表示IES1离网
t_realtime = zeros(24*period,3);
isDAsingle=0;
if isDAsingle==0
    subgrad_dayahead_180729;
end
%用次梯度法解日前优化问题

for t_current = 1:24*period
    disp(['t_current is ',num2str(t_current)]);
    if isDAsingle==1
   
    if t_current == 1 %则进行日前优化
        tic
        % 负荷和可再生能源的预测值没有变化
%         EH1.predict(0);
%         EH2.predict(0);
%         EH3.predict(0);
        for pt = t_current : 1 : 24*period
            tic
%             Grid1.zGenerate(elePrice(pt)); % 在当前优化周期内不变

            % 最低价格，一般都是需求大于供给
            priceArray(pt) = minMarketPrice;
            [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, pt);
            clearDemand_minPrice_EH1 = x(1);
            [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, pt);
            clearDemand_minPrice_EH2 = x(1);
            [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice3, pt);
            clearDemand_minPrice_EH3 = x(1);   
            clearDemand_minPrice_grid = Grid1.handlePrice(priceArray(pt), pt);
            clearDemand_minPrice = [clearDemand_minPrice_grid; clearDemand_minPrice_EH1; clearDemand_minPrice_EH2; clearDemand_minPrice_EH3]; % 需求为正，供给为负
            
            % 最高价格，一般都是供给大于需求
            priceArray(pt) = maxMarketPrice;
            [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, pt);
            clearDemand_maxPrice_EH1 = x(1);
            [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, pt);
            clearDemand_maxPrice_EH2 = x(1);
            [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice3, pt);
            clearDemand_maxPrice_EH3 = x(1);   
            clearDemand_maxPrice_grid = Grid1.handlePrice(priceArray(pt), pt);
            clearDemand_maxPrice = [clearDemand_maxPrice_grid; clearDemand_maxPrice_EH1; clearDemand_maxPrice_EH2; clearDemand_maxPrice_EH3]; % 需求为正，供给为负
            
            iterationNumber = 2;
            if sum(clearDemand_minPrice) * sum(clearDemand_maxPrice) <= 0 % 说明出清点在这个区间内，有两个问题，一是等于零是否直接结束，二是如果出清点不唯一怎么办
                % 市场出清得到出清价格，并更新预测电价序列
                [priceArray(pt), clearDemand] = iterativeClear(minMarketPrice, maxMarketPrice, clearDemand_minPrice, clearDemand_maxPrice, ee, priceArray, gasPrice1, gasPrice3, pt);
            else
                disp('Clearing point is not in the given interval.')
            end
            iterationTimes(pt,1) = iterationNumber;
            
            % 得到出清价格后，还要明确出清功率（EH更新自身状态），但此时不能保证各出清功率之和为零？
            gridClearDemand(pt) = clearDemand(1);
            EH1.conditionHandlePrice_2(priceArray, gasPrice1, pt, clearDemand(2));
            EH2.conditionHandlePrice_2(priceArray, gasPrice1, pt, clearDemand(3));
            EH3.conditionHandlePrice_2(priceArray, gasPrice3, pt, clearDemand(4));
            
        end
        
        priceArray_record(:,1) = priceArray;
        t_dayahead = toc; %这个时间不准确，因为3个IES应该是并行计算的
        
        % 日前优化的结果
        [result_Ele(:,1), result_CHP_G(:,1), result_Boiler_G(:,1), result_ES_discharge(:,1), result_ES_charge(:,1), result_HS_discharge(:,1), result_HS_charge(:,1), result_ES_SOC(:,1), result_HS_SOC(:,1), EH1_Le, EH1_Lh, EH1_solarP, EH1_windP] = EH1.getResult;
        [result_Ele(:,2), result_CHP_G(:,2), result_Boiler_G(:,2), result_ES_discharge(:,2), result_ES_charge(:,2), result_HS_discharge(:,2), result_HS_charge(:,2), result_ES_SOC(:,2), result_HS_SOC(:,2), EH2_Le, EH2_Lh, EH2_solarP, EH2_windP] = EH2.getResult;
        [result_Ele(:,3), result_CHP_G(:,3), result_Boiler_G(:,3), result_ES_discharge(:,3), result_ES_charge(:,3), result_HS_discharge(:,3), result_HS_charge(:,3), result_ES_SOC(:,3), result_HS_SOC(:,3), EH3_Le, EH3_Lh, ~, EH3_windP] = EH3.getResult;
    end
    clearDemand_grid_sin=sum(result_Ele');    
    priceArray_pre_sin=priceArray;
    end

    off_grid = 0;
    
    % 日内优化，分为正常运行和EH1离网   
    if t_current==11
       a=1; 
    end
    if off_grid == 0 % 正常运行
        EH1.predict(t_current);
        EH2.predict(t_current);
        EH3.predict(t_current);
        
%         Grid1.zGenerate(elePrice(t_current)); %在当前优化周期内不变
        
        % 最低价格，一般都是需求大于供给
        priceArray(t_current) = minMarketPrice;
        [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, t_current);
        clearDemand_minPrice_EH1 = x(1);
        [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, t_current);
        clearDemand_minPrice_EH2 = x(1);
        [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice3, t_current);
        clearDemand_minPrice_EH3 = x(1);
        clearDemand_minPrice_grid = Grid1.handlePrice(priceArray(t_current), t_current);
        clearDemand_minPrice = [clearDemand_minPrice_grid; clearDemand_minPrice_EH1; clearDemand_minPrice_EH2; clearDemand_minPrice_EH3]; % 需求为正，供给为负
            
        % 最高价格，一般都是供给大于需求
        priceArray(t_current) = maxMarketPrice;
        [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, t_current);
        clearDemand_maxPrice_EH1 = x(1);
        [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, t_current);
        clearDemand_maxPrice_EH2 = x(1);
        [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice3, t_current);
        clearDemand_maxPrice_EH3 = x(1);
        clearDemand_maxPrice_grid = Grid1.handlePrice(priceArray(t_current), t_current);
        clearDemand_maxPrice = [clearDemand_maxPrice_grid; clearDemand_maxPrice_EH1; clearDemand_maxPrice_EH2; clearDemand_maxPrice_EH3]; % 需求为正，供给为负
            
        iterationNumber = 2;
        if sum(clearDemand_minPrice) * sum(clearDemand_maxPrice) <= 0 % 说明出清点在这个区间内，有两个问题，一是等于零是否直接结束，二是如果出清点不唯一怎么办
            %市场出清得到出清价格，并更新预测电价序列
            [priceArray(t_current), clearDemand] = iterativeClear(minMarketPrice, maxMarketPrice, clearDemand_minPrice, clearDemand_maxPrice, ee, priceArray, gasPrice1, gasPrice3, t_current);
        else
            disp('Clearing point is not in the given interval.')
        end
        iterationTimes(t_current,2) = iterationNumber;
        
        % 得到出清价格后，还要明确出清功率（EH更新自身状态）
        gridClearDemand(t_current) = clearDemand(1);
        EH1.conditionHandlePrice_2(priceArray, gasPrice1, t_current, clearDemand(2));
        EH2.conditionHandlePrice_2(priceArray, gasPrice1, t_current, clearDemand(3));
        EH3.conditionHandlePrice_2(priceArray, gasPrice3, t_current, clearDemand(4));
        
    else % IES1离网
        
%         EH1.predict(t_current);
        EH2.predict(t_current);
        EH3.predict(t_current);
        
%         Grid1.zGenerate(elePrice(t_current)); %在当前优化周期内不变
        
        % 最低价格，一般都是需求大于供给
        priceArray(t_current) = minMarketPrice;
%         [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, t_current);
%         clearDemand_minPrice_EH1 = x(1);
        [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, t_current);
        clearDemand_minPrice_EH2 = x(1);
        [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice3, t_current);
        clearDemand_minPrice_EH3 = x(1);
        clearDemand_minPrice_grid = Grid1.handlePrice(priceArray(t_current), t_current);
        clearDemand_minPrice = [clearDemand_minPrice_grid; clearDemand_minPrice_EH2; clearDemand_minPrice_EH3]; % 需求为正，供给为负
            
        % 最高价格，一般都是供给大于需求
        priceArray(t_current) = maxMarketPrice;
%         [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, t_current);
%         clearDemand_maxPrice_EH1 = x(1);
        [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, t_current);
        clearDemand_maxPrice_EH2 = x(1);
        [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice3, t_current);
        clearDemand_maxPrice_EH3 = x(1);
        clearDemand_maxPrice_grid = Grid1.handlePrice(priceArray(t_current), t_current);
        clearDemand_maxPrice = [clearDemand_maxPrice_grid; clearDemand_maxPrice_EH2; clearDemand_maxPrice_EH3]; % 需求为正，供给为负
            
        iterationNumber = 2;
        if sum(clearDemand_minPrice) * sum(clearDemand_maxPrice) <= 0 % 说明出清点在这个区间内，有两个问题，一是等于零是否直接结束，二是如果出清点不唯一怎么办
            %市场出清得到出清价格，并更新预测电价序列
            [priceArray(t_current), clearDemand] = iterativeClear(minMarketPrice, maxMarketPrice, clearDemand_minPrice, clearDemand_maxPrice, ee, priceArray, gasPrice1, gasPrice3, t_current);
        else
            disp('Clearing point is not in the given interval.')
        end
        iterationTimes(t_current,2) = iterationNumber;
        
        % 得到出清价格后，还要明确出清功率（EH更新自身状态）
        gridClearDemand(t_current) = clearDemand(1);
%         EH1.conditionHandlePrice_2(priceArray, gasPrice1, t_current, clearDemand(2));
        EH2.conditionHandlePrice_2(priceArray, gasPrice1, t_current, clearDemand(2));
        EH3.conditionHandlePrice_2(priceArray, gasPrice3, t_current, clearDemand(3));
    
    end
end
priceArray_record(:,2) = priceArray;

[result_Ele(:,1), result_CHP_G(:,1), result_Boiler_G(:,1), result_ES_discharge(:,1), result_ES_charge(:,1), result_HS_discharge(:,1), result_HS_charge(:,1), result_ES_SOC(:,1), result_HS_SOC(:,1), EH1_Le, EH1_Lh, EH1_solarP, EH1_windP] = EH1.getResult;
[result_Ele(:,2), result_CHP_G(:,2), result_Boiler_G(:,2), result_ES_discharge(:,2), result_ES_charge(:,2), result_HS_discharge(:,2), result_HS_charge(:,2), result_ES_SOC(:,2), result_HS_SOC(:,2), EH2_Le, EH2_Lh, EH2_solarP, EH2_windP] = EH2.getResult;
[result_Ele(:,3), result_CHP_G(:,3), result_Boiler_G(:,3), result_ES_discharge(:,3), result_ES_charge(:,3), result_HS_discharge(:,3), result_HS_charge(:,3), result_ES_SOC(:,3), result_HS_SOC(:,3), EH3_Le, EH3_Lh, EH3_solarP, EH3_windP] = EH3.getResult;



%迭代出清，梯度法
%{
priceArray = elePrice; %由历史数据得到预测电价 %也是用于迭代的价格变量
priceArray_record = zeros(24*period,2); %一列日前，一列实时
ee = 0.01; %0.0001 0.0003
iterativeStep = 0.000001; %0.00001 0.0001
iterationTimes = zeros(24*period, 2); %记录迭代次数
maxIteration = 3000; %最大迭代次数

demand_sum = zeros(priceNumbers, 1);
gridClearDemand = zeros(24*period,1);

off_grid = 0; % 0表示正常运行，1表示IES1离网
t_realtime = zeros(24*period,3);

% 用于测试的电价
% priceArray = [0.334979896750537;0.345365343857041;0.332393081385998;0.308531501102345;0.322635963415099;0.304882110539400;0.304515835848279;0.304153787790538;1.04508174057890;1.04470000000178;1.04470000000007;1.04470000000001;0.626800000000104;0.626800000000167;0.626800000055620;0.627101953428408;0.628016602128634;1.04470000000009;1.04512412592951;1.04574079511051;1.04623511011097;0.627798249878306;0.641181118337553;0.633909463862792];


for t_current = 1:24*period
    disp(['t_current is ',num2str(t_current)]);
    
    if t_current == 1 %则进行日前优化
        tic
%         EH1.predict(0);
%         EH2.predict(0);
%         EH3.predict(0);
        for pt = t_current : 1 : 24*period
            tic
            
%             figure(1);hold on;

            number = 1;
            
            lamda_old = -10;
            lamda_new = 0.0; %取初始值：对预测电价没有偏差
            lamda_record = zeros(maxIteration+1, 1);            
            lamda_record(number) = lamda_new;
            
            lamda_avg_old = lamda_old;
            lamda_avg_new = lamda_new;
            lamda_avg_record = zeros(maxIteration+1, 1);
            lamda_avg_record(number) = lamda_avg_new;
            
            clearDemand_record = zeros(maxIteration+1, 1);
        
            
            %如果前后两次价格的偏差太大，则返回第1步
            while number<=2 || abs(lamda_avg_new - lamda_avg_old) > ee || sum(clearDemand_new) * sum(clearDemand_old) > 1e-4  %1e-6, 不能直接取0
                % 后一个条件是因为即使lamda收敛后，供需也不平衡，所以需要取一正一负两个点，来求零点
                % && || 的前一个为否，则后一个就不计算了
                % 要求至少迭代两次（number=1，2）
                
                if number > maxIteration
                    error('超出最大迭代次数');
                end
                
                %当前价格下的出力
                priceArray(pt) = elePrice(pt) + lamda_new;
                [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, pt);
                clearDemand_EH1_new = x(1);
                [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, pt);
                clearDemand_EH2_new = x(1);
                [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice3, pt);
                clearDemand_EH3_new = x(1);
                
%                 clearDemand_grid_new = Grid1.handlePrice(priceArray(pt), pt);

                f1 = - lamda_new;
                lb1 = eleLimit_total(2);
                ub1 = eleLimit_total(1);
                [clearDemand_grid_new, value1, flag1]  = linprog(f1, [], [], [], [], lb1, ub1);    
                
                % 存储老的clearDemand，计算新的clearDemand，并记录
                if number>1
                    clearDemand_old = clearDemand_new; % number=2时才记录第一次
                end
                clearDemand_new = [-clearDemand_grid_new; clearDemand_EH1_new; clearDemand_EH2_new; clearDemand_EH3_new]; % 需求为正，供给为负
                clearDemand_record(number) = sum(clearDemand_new);
                
                % 存储老的lamda，计算新的lamda（通过梯度法），并记录
                lamda_old = lamda_new;
%                 lamda_new = max(0, lamda_old + sum(clearDemand) * iterativeStep);
                lamda_new = lamda_old + sum(clearDemand_new) * iterativeStep;
                number = number + 1;
                lamda_record(number) = lamda_new;

                % 存储老的平均值，计算新的平均值，并记录
                lamda_avg_old = lamda_avg_new;
                % 系数bb：新值占平均值的比例
                bb = 1/number; % 方法一：平均值，实际就是 = (lamda_avg_old * (number-1) + lamda_new * 1) / number;
                % bb = 0.5; % 方法二：加权平均，重视新值，试了结果不能收敛
                lamda_avg_new = lamda_avg_old * (1-bb) + lamda_new * bb; 
                lamda_avg_record(number) = lamda_avg_new;
                
%                 lamdaArrayAvg_new = sum(lamdaArray) / length(lamdaArray); %方法一
%                 lamda_avg_new = 1 / length(lamda_record) * lamda_record(number) + (length(lamda_record)-1) / length(lamda_record) * lamda_avg_old; %方法二                    
            end
            
            % 否则停止迭代
            iterationTimes(pt,1) = number - 1;
            % 得到出清电价，假设线性，根据迭代最后两次的结果，求新的出清价格和出清功率
            if sum(clearDemand_new) * sum(clearDemand_old) <= 0
                if lamda_record(number-1) == lamda_record(number-2) %防止最后的计算式的分母为零
                    clearLamda = lamda_record(number-1);
                elseif clearDemand_record(number-1) == clearDemand_record(number-2) %防止最后的计算式的分母为零
                    clearLamda = (lamda_record(number-1) + lamda_record(number-2)) / 2;
                else
                    slope = (clearDemand_record(number-1) - clearDemand_record(number-2)) / (lamda_record(number-1) - lamda_record(number-2));
                    clearLamda = lamda_record(number-2) + (0 - clearDemand_record(number-2)) / slope;
                end
            else % 值在[0，1e-4]之间，那么就没有零点了
                clearLamda = (lamda_record(number-1) + lamda_record(number-2)) / 2;
            end
            
            clearDemand = zeros(length(clearDemand_new) ,1);
            for i=1:length(clearDemand_new)
                if lamda_record(number-1) == lamda_record(number-2)
                    clearDemand(i) = clearDemand_new(i);
                else
                    slope = (clearDemand_new(i) - clearDemand_old(i)) / (lamda_record(number-1) - lamda_record(number-2));
                    clearDemand(i) = clearDemand_new(i) + (clearLamda - lamda_record(number-1)) * slope;
                end
            end
            
            % 现在出清电价可正可负
%             if clearPrice < minMarketPrice || clearPrice > maxMarketPrice
%                 error('出清电价超出允许范围')
%             end
                       
            % 根据得到的出清价格以及出清功率，EH进行一次优化，以更新自身状态
            priceArray(pt) = elePrice(pt) + clearLamda;
            gridClearDemand(pt) = clearDemand(1);
            EH1.conditionHandlePrice_2(priceArray, gasPrice1, pt, clearDemand(2));
            EH2.conditionHandlePrice_2(priceArray, gasPrice1, pt, clearDemand(3));
            EH3.conditionHandlePrice_2(priceArray, gasPrice3, pt, clearDemand(4));
            
        end
        
        priceArray_record(:,1) = priceArray;
        t_dayahead = toc; %这个时间不准确，因为3个IES应该是并行计算的
        
        % 日前优化的结果
        [result_Ele(:,1), result_CHP_G(:,1), result_Boiler_G(:,1), result_ES_discharge(:,1), result_ES_charge(:,1), result_HS_discharge(:,1), result_HS_charge(:,1), result_ES_SOC(:,1), result_HS_SOC(:,1), EH1_Le, EH1_Lh, EH1_solarP, EH1_windP] = EH1.getResult;
        [result_Ele(:,2), result_CHP_G(:,2), result_Boiler_G(:,2), result_ES_discharge(:,2), result_ES_charge(:,2), result_HS_discharge(:,2), result_HS_charge(:,2), result_ES_SOC(:,2), result_HS_SOC(:,2), EH2_Le, EH2_Lh, EH2_solarP, EH2_windP] = EH2.getResult;
        [result_Ele(:,3), result_CHP_G(:,3), result_Boiler_G(:,3), result_ES_discharge(:,3), result_ES_charge(:,3), result_HS_discharge(:,3), result_HS_charge(:,3), result_ES_SOC(:,3), result_HS_SOC(:,3), EH3_Le, EH3_Lh, EH3_solarP, EH3_windP] = EH3.getResult;
    end
    
    
    
    % 日内优化，分为正常运行和EH1离网   
    if off_grid == 0 % 正常运行
        EH1.predict(t_current);
        EH2.predict(t_current);
        EH3.predict(t_current);
        
        number = 1;
        
        lamda_old = -10;
        lamda_new = priceArray_record(t_current,1) - elePrice(t_current); %取初始值：为日前优化的出清lamda，当时没保存，需要计算一下
        lamda_record = zeros(maxIteration+1, 1);
        lamda_record(number) = lamda_new;
        
        lamda_avg_old = lamda_old;
        lamda_avg_new = lamda_new;
        lamda_avg_record = zeros(maxIteration+1, 1);
        lamda_avg_record(number) = lamda_avg_new;
        
        clearDemand_record = zeros(maxIteration+1, 1);
        
        %如果前后两次价格的偏差太大，则返回第1步
        while number<=2 || abs(lamda_avg_new - lamda_avg_old) > ee || sum(clearDemand_new) * sum(clearDemand_old) > 1e-4  %1e-6, 不能直接取0
            % 后一个条件是因为即使lamda收敛后，供需也不平衡，所以需要取一正一负两个点，来求零点
            % && || 的前一个为否，则后一个就不计算了
            % 要求至少迭代两次（number=1，2）     
            
            if number > maxIteration
                error('超出最大迭代次数');
            end
            
            %当前价格下的出力
            priceArray(t_current) = elePrice(t_current) + lamda_new;
            [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, t_current);
            clearDemand_EH1_new = x(1);
            [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, t_current);
            clearDemand_EH2_new = x(1);
            [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice3, t_current);
            clearDemand_EH3_new = x(1);
            
            f1 = - lamda_new;
            lb1 = eleLimit_total(2);
            ub1 = eleLimit_total(1);
            [clearDemand_grid_new, value1, flag1]  = linprog(f1, [], [], [], [], lb1, ub1);
            
            % 存储老的clearDemand，计算新的clearDemand，并记录
            if number>1
                clearDemand_old = clearDemand_new;
            end
            clearDemand_new = [-clearDemand_grid_new; clearDemand_EH1_new; clearDemand_EH2_new; clearDemand_EH3_new]; % 需求为正，供给为负
            clearDemand_record(number) = sum(clearDemand_new);
            
            % 存储老的lamda，计算新的lamda（通过梯度法），并记录
            lamda_old = lamda_new;
            lamda_new = lamda_old + sum(clearDemand_new) * iterativeStep;
            number = number + 1;
            lamda_record(number) = lamda_new;
            
            % 存储老的平均值，计算新的平均值，并记录
            lamda_avg_old = lamda_avg_new;
            % 系数bb：新值占平均值的比例
            bb = 1/number; % 方法一：平均值，实际就是 = (lamda_avg_old * (number-1) + lamda_new * 1) / number;
            % bb = 0.5; % 方法二：加权平均，重视新值，试了结果不能收敛
            lamda_avg_new = lamda_avg_old * (1-bb) + lamda_new * bb;
            lamda_avg_record(number) = lamda_avg_new;
        end
        
        % 否则停止迭代
        iterationTimes(t_current,2) = number - 1;
        % 得到出清电价，假设线性，根据迭代最后两次的结果，求新的出清价格和出清功率
        if sum(clearDemand_new) * sum(clearDemand_old) <= 0
            if lamda_record(number-1) == lamda_record(number-2) %防止最后的计算式的分母为零
                clearLamda = lamda_record(number-1);
            elseif clearDemand_record(number-1) == clearDemand_record(number-2) %防止最后的计算式的分母为零
                clearLamda = (lamda_record(number-1) + lamda_record(number-2)) / 2;
            else
                slope = (clearDemand_record(number-1) - clearDemand_record(number-2)) / (lamda_record(number-1) - lamda_record(number-2));
                clearLamda = lamda_record(number-2) + (0 - clearDemand_record(number-2)) / slope;
            end
        else % 值在[0，1e-4]之间，那么就没有零点了
            clearLamda = (lamda_record(number-1) + lamda_record(number-2)) / 2;
        end
        
        clearDemand = zeros(length(clearDemand_new) ,1);
        for i=1:length(clearDemand_new)
            if lamda_record(number-1) == lamda_record(number-2)
                clearDemand(i) = clearDemand_new(i);
            else
                slope = (clearDemand_new(i) - clearDemand_old(i)) / (lamda_record(number-1) - lamda_record(number-2));
                clearDemand(i) = clearDemand_new(i) + (clearLamda - lamda_record(number-1)) * slope;
            end
        end
        
        % 根据得到的出清价格以及出清功率，EH进行一次优化，以更新自身状态
        priceArray(t_current) = elePrice(t_current) + clearLamda;
        gridClearDemand(t_current) = clearDemand(1);
        EH1.conditionHandlePrice_2(priceArray, gasPrice1, t_current, clearDemand(2));
        EH2.conditionHandlePrice_2(priceArray, gasPrice1, t_current, clearDemand(3));
        EH3.conditionHandlePrice_2(priceArray, gasPrice3, t_current, clearDemand(4)); 
        
    else % IES1离网
    end
end
priceArray_record(:,2) = priceArray;

[result_Ele(:,1), result_CHP_G(:,1), result_Boiler_G(:,1), result_ES_discharge(:,1), result_ES_charge(:,1), result_HS_discharge(:,1), result_HS_charge(:,1), result_ES_SOC(:,1), result_HS_SOC(:,1), EH1_Le, EH1_Lh, EH1_solarP, EH1_windP] = EH1.getResult;
[result_Ele(:,2), result_CHP_G(:,2), result_Boiler_G(:,2), result_ES_discharge(:,2), result_ES_charge(:,2), result_HS_discharge(:,2), result_HS_charge(:,2), result_ES_SOC(:,2), result_HS_SOC(:,2), EH2_Le, EH2_Lh, EH2_solarP, EH2_windP] = EH2.getResult;
[result_Ele(:,3), result_CHP_G(:,3), result_Boiler_G(:,3), result_ES_discharge(:,3), result_ES_charge(:,3), result_HS_discharge(:,3), result_HS_charge(:,3), result_ES_SOC(:,3), result_HS_SOC(:,3), EH3_Le, EH3_Lh, EH3_solarP, EH3_windP] = EH3.getResult;
%}



