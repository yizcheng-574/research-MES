% 集中式优化：单次
clc;clear;
caseType = 2;
para_init;
%1 全时段集中式优化
%2 集中式滚动优化
%0 分布式滚动优化
global period
for pt = 1: 24 * period
    EH1.predict(pt);
    EH2.predict(pt);
    EH3.predict(pt);
    [EH1_f, EH1_intcon, EH1_A, EH1_b, EH1_Aeq, EH1_beq, EH1_lb, EH1_ub, EH1_A_eleLimit_total] = EH1.MilpMatrix(gasPrice1, pt);
    [EH2_f, EH2_intcon, EH2_A, EH2_b, EH2_Aeq, EH2_beq, EH2_lb, EH2_ub, EH2_A_eleLimit_total] = EH2.MilpMatrix(gasPrice1, pt);
    [EH3_f, EH3_intcon, EH3_A, EH3_b, EH3_Aeq, EH3_beq, EH3_lb, EH3_ub, EH3_A_eleLimit_total] = EH3.MilpMatrix(gasPrice1, pt);
    
    time = 24 * period - pt + 1; %总时间段
    number = IESNUMBER;
    var = time * 13;
    totalVar = time * 13 * number; %总变量数
    %第1,2,3组time是购电量、CHP购气量、锅炉购气量，第4-7组time是储电、储热的放、充功率
    %第8-9组是可平移负荷，第10-11组是EESbinary， 第12-13组是TESbinary
    intcon = [EH1_intcon, EH2_intcon + var, EH3_intcon + var];
    f = [EH1_f; EH2_f; EH3_f];
    ub = [EH1_ub; EH2_ub; EH3_ub];
    lb = [EH1_lb; EH2_lb; EH3_lb];
    [lr , lc] = size(EH1_Aeq);
    Aeq = zeros( lr*number, lc*number);
    for i= 1 : number
        eval(['Aeq((i-1)*lr+1:lr*i,(i-1)*lc+1:lc*i) = EH',num2str(i),'_Aeq;']);
    end
    beq = [EH1_beq; EH2_beq; EH3_beq];
    
    [lr , lc] = size(EH1_A);
    A1 = zeros( lr*number, lc*number);
    for i= 1 : number
        eval(['A1((i-1)*lr+1:lr*i,(i-1)*lc+1:lc*i) = EH',num2str(i),'_A;']);
    end
    b1 = [EH1_b; EH2_b;EH3_b];
    
    % 需要额外增加一个购电量的上、下限约束
    A2 = [EH1_A_eleLimit_total, EH2_A_eleLimit_total, EH2_A_eleLimit_total];
    b2 = ones(time, 1) .* eleLimit_total(1);
%     b2_sale = ones(time, 1) .* eleLimit_total(2);
    b2_sale = zeros(time, 1);
    A = [A1; A2; -A2];
    b = [b1; b2; -b2_sale];
    
    [x,fval,exitflag,output] = intlinprog(f,intcon, A, b, Aeq, beq, lb, ub);
    
    for IES_no = 1 : IESNUMBER
        eval(['EH',num2str(IES_no),'.update_central(x, pt, IES_no, 1);']);
    end
    
end

[result_Ele(:,4), result_CHP_G(:,4), result_Boiler_G(:,4), result_ES_discharge(:,4), result_ES_charge(:,4), result_HS_discharge(:,4), result_HS_charge(:,4), result_ES_SOC(:,4), result_HS_SOC(:,4), result_EH_Le(:, 4), result_EH_Lh(:,4), result_EH_solarP(:,4), result_EH_windP(:,4), result_EH_Edr(:,4), result_EH_Hdr(:,4)] = EH1.getResult;
[result_Ele(:,5), result_CHP_G(:,5), result_Boiler_G(:,5), result_ES_discharge(:,5), result_ES_charge(:,5), result_HS_discharge(:,5), result_HS_charge(:,5), result_ES_SOC(:,5), result_HS_SOC(:,5), result_EH_Le(:, 5), result_EH_Lh(:,5), result_EH_solarP(:,5), result_EH_windP(:,5), result_EH_Edr(:,5), result_EH_Hdr(:,5)] = EH2.getResult;
[result_Ele(:,6), result_CHP_G(:,6), result_Boiler_G(:,6), result_ES_discharge(:,6), result_ES_charge(:,6), result_HS_discharge(:,6), result_HS_charge(:,6), result_ES_SOC(:,6), result_HS_SOC(:,6), result_EH_Le(:, 6), result_EH_Lh(:,6), result_EH_solarP(:,6), result_EH_windP(:,6), result_EH_Edr(:,6), result_EH_Hdr(:,6)] = EH3.getResult;

main_central;
main_handle_MILP;