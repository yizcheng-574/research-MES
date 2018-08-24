% clc
% clear

%多能协调优化
%能源枢纽的本地集中式优化策略
%根据输入参数，形成所需要的矩阵

%20171010 加可再生能源，改为滚动优化
%20180110 储能平衡性约束改为不等式约束，这样就没有等式约束了
%20180129 增加购电效率，即线路传输的线损率Ele_eff，增加period

function [f, ub, lb, A, b, A_eleLimit_total] = OptMatrix_rolling_20171010(eleLimit, gasLimit, Le, Lh, CHP_para, Boiler_para, ES_para, HS_para, Gprice, windP, solarP, t_current)
global period elePrice

%联络线电交换和气交换约束
Ele_max = eleLimit(1);
Ele_min = eleLimit(2);
Ele_eff = eleLimit(3); %购电效率，即线路传输的线损率
Gas_max = gasLimit; %两个购气量的和
% Gas_min = 0;

%CHP的参数
CHP_GE_eff = CHP_para(1);
CHP_GH_eff = CHP_para(2);
CHP_G_max = CHP_para(3) / CHP_GE_eff; %由额定电功率得到最大购气量
CHP_G_min = CHP_para(4) * CHP_G_max;

%锅炉
Boiler_eff = Boiler_para(1);
Boiler_G_max = Boiler_para(2) / Boiler_eff;  %10
% Boiler_G_min = 0;

%电储能和热储能
ES_totalC = ES_para(1);
ES_maxSOC = ES_para(2);
ES_minSOC = ES_para(3);
ES_currentSOC = ES_para(4);
ES_targetSOC = ES_para(5);
ES_Pmax = ES_totalC / ES_para(6);
ES_eff = ES_para(7);

HS_totalC = HS_para(1);
HS_maxSOC = HS_para(2);
HS_minSOC = HS_para(3);
HS_currentSOC = HS_para(4);
HS_targetSOC = HS_para(5);
HS_Hmax = HS_totalC / HS_para(6);
HS_eff = HS_para(7);



time = 24*period - t_current + 1; %总时间段
var = time * 7; %总变量数
%第1,2,3组time是购电量、CHP购气量、锅炉购气量，第4-7组time是储电、储热的放、充功率



%求一次项系数f 是个列向量
f = zeros(var, 1);
for i = 1 : time
    % f(i, 1) = elePrice(i);
    f(i, 1) = elePrice(t_current + i - 1); %elePrice的size不变，但是取部分值
    f(time+i, 1) = Gprice;
    f(time*2+i, 1) = Gprice;
end

%变量上下限
ub = zeros(var, 1);
lb = zeros(var, 1);
for i = 1 : time
    ub(i, 1) = Ele_max;
    ub(time+i, 1) = CHP_G_max;
    ub(time*2+i, 1) = Boiler_G_max;
    ub(time*3+i, 1) = ES_Pmax;
    ub(time*4+i, 1) = ES_Pmax;
    ub(time*5+i, 1) = HS_Hmax;
    ub(time*6+i, 1) = HS_Hmax;
end
for i = 1 : time
    lb(i, 1) = Ele_min;
    lb(time+i, 1) = CHP_G_min;
    %     lb(time*2+i, 1) = 0;
    %     lb(time*3+i, 1) = 0;
    %     lb(time*4+i, 1) = 0;
    %     lb(time*5+i, 1) = 0;
    %     lb(time*6+i, 1) = 0;
end



%等式约束包括：无


%不等式约束包括：储能平衡性约束；功率平衡约束（供大于求）；SOC约束；购气量和的约束；（爬坡率约束）
%电、热储能平衡性约束（改为不等式约束 截止时刻值>目标值 即 -截止时刻值<-目标值）
Aeq_ES = zeros(1, var);
beq_ES = - (ES_targetSOC - ES_currentSOC) * ES_totalC;
for i=1:time
    Aeq_ES(1, time*3+i) = 1/ES_eff; %放电
    Aeq_ES(1, time*4+i) = - 1*ES_eff; %充电
end
Aeq_HS = zeros(1, var);
beq_HS = - (HS_targetSOC - HS_currentSOC) * HS_totalC;
for i=1:time
    Aeq_HS(1, time*5+i) = 1/HS_eff; %放热
    Aeq_HS(1, time*6+i) = - 1*HS_eff; %充热
end

%电、热平衡约束
Aeq_Ebus = zeros(time, var);
Aeq_Hbus = zeros(time, var);
% beq_Ebus = - Le;
% beq_Hbus = - Lh;
beq_Ebus = - Le(t_current : 24*period) + windP(t_current : 24*period) + solarP(t_current : 24*period); %Le的size不变，但是取部分值
beq_Hbus = - Lh(t_current : 24*period); %Lh的size不变，但是取部分值 

for i=1:time
    Aeq_Ebus(i,i) = - Ele_eff;
    Aeq_Ebus(i,time+i) = - CHP_GE_eff;
    Aeq_Ebus(i,time*3+i) = - 1; %放电
    Aeq_Ebus(i,time*4+i) = 1; %充电
end
for i=1:time
    Aeq_Hbus(i,time+i) = - CHP_GH_eff;
    Aeq_Hbus(i,time*2+i) = - Boiler_eff;
    Aeq_Hbus(i,time*5+i) = - 1; %放热
    Aeq_Hbus(i,time*6+i) = 1; %充热
end

%SOC约束 A1是上限，A2是下限
A1_Esoc = zeros(time, var);
A2_Esoc = zeros(time, var);
b1_Esoc = ones(time,1) * (ES_maxSOC - ES_currentSOC) * ES_totalC;
b2_Esoc = ones(time,1) * (ES_currentSOC - ES_minSOC) * ES_totalC;
for i=1:time
    for j=1 : i
        A1_Esoc(i, time*3+j) = -1/ES_eff; %放电
        A1_Esoc(i, time*4+j) = 1*ES_eff; %充电
    end
end
for i=1:time
    for j=1 : i
        A2_Esoc(i, time*3+j) = 1/ES_eff; %放电
        A2_Esoc(i, time*4+j) = -1*ES_eff; %充电
    end
end

A1_Hsoc = zeros(time, var);
A2_Hsoc = zeros(time, var);
b1_Hsoc = ones(time,1) * (HS_maxSOC - HS_currentSOC) * HS_totalC;
b2_Hsoc = ones(time,1) * (HS_currentSOC - HS_minSOC) * HS_totalC;
for i=1:time
    for j=1 : i
        A1_Hsoc(i, time*5+j) = -1/HS_eff; %放热
        A1_Hsoc(i, time*6+j) = 1*HS_eff; %充热
    end
end
for i=1:time
    for j=1 : i
        A2_Hsoc(i, time*5+j) = 1/HS_eff; %放热
        A2_Hsoc(i, time*6+j) = -1*HS_eff; %充热
    end
end

%购气量和的约束
A_Gmax = zeros(time, var);
b_Gmax = ones(time,1) .* Gas_max;
for i=1:time
    A_Gmax(i, time+i) = 1;
    A_Gmax(i, time*2+i) = 1;
end



%归纳所有线性约束
%不等式约束包括：储能平衡性约束，功率平衡约束，SOC约束，购气量和的约束，（爬坡率约束）
A=[Aeq_Ebus; Aeq_Hbus;   Aeq_ES; Aeq_HS;   A1_Esoc; A2_Esoc; A1_Hsoc; A2_Hsoc; A_Gmax];
b=[beq_Ebus; beq_Hbus;   beq_ES; beq_HS;   b1_Esoc; b2_Esoc; b1_Hsoc; b2_Hsoc; b_Gmax];


% 需要额外增加一个购电量的上、下限约束
A_eleLimit_total = zeros(time, var);
for i=1:time
    A_eleLimit_total(i, i) = 1;
end

end




