% 集中式优化：单次
clc;clear;
para_init;
isCentral = 1;
%1 全时段集中式优化
%2 集中式滚动优化
%0 分布式滚动优化
global period
if isCentral == 2 
    temporal = 24* period;
else
    temporal = 1;
%     EH1.predict(0);
%     EH2.predict(0);
%     EH3.predict(0);
end
for pt = 1: temporal
    EH1.predict(pt);
    EH2.predict(pt);
    EH3.predict(pt);
    [EH1_f, EH1_ub, EH1_lb, EH1_Aeq, EH1_beq, EH1_A, EH1_b, EH1_A_eleLimit_total ] = EH1.OptMatrix(gasPrice1, pt);
    [EH2_f, EH2_ub, EH2_lb, EH2_Aeq, EH2_beq, EH2_A, EH2_b, EH2_A_eleLimit_total ] = EH2.OptMatrix(gasPrice1, pt);
    [EH3_f, EH3_ub, EH3_lb, EH3_Aeq, EH3_beq, EH3_A, EH3_b, EH3_A_eleLimit_total ] = EH3.OptMatrix(gasPrice1, pt);
    time = 24*period - pt + 1; %总时间段
    number = IESNUMBER;
    var = time * 9;
    totalVar = time * 9 * number; %总变量数
    %第1,2,3组time是购电量、CHP购气量、锅炉购气量，第4-7组time是储电、储热的放、充功率
    
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
    b2_sale = ones(time, 1) .* eleLimit_total(2);
    A = [A1; A2; -A2];
    b = [b1; b2; -b2_sale];
    
    
    [x,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub);
    if isCentral == 2
        for IES_no = 1 : IESNUMBER
            eval(['EH',num2str(IES_no),'.update_central(x, pt, IES_no);']);
        end
    end
end
if isCentral == 1 
    for IES_no = 1 : number
        result_Ele(:, IES_no) = x(1 + (IES_no - 1)* var :time  + (IES_no - 1)* var);
        result_CHP_G(:,IES_no) = x(time + 1 + (IES_no - 1)* var :time * 2 +  (IES_no - 1)* var);
        result_Boiler_G(:,IES_no) = x(time * 2 + 1 + (IES_no - 1)* var :time * 3 +  (IES_no - 1)* var);
        result_ES_discharge(:,IES_no) = x(time * 3 + 1 + (IES_no - 1)* var :time * 4 +  (IES_no - 1)* var);
        result_ES_charge(:,IES_no) = x(time * 4 + 1 + (IES_no - 1)* var :time * 5 +  (IES_no - 1)* var);
        result_HS_discharge(:,IES_no) = x(time * 5 + 1 + (IES_no - 1)* var :time * 6 +  (IES_no - 1)* var);
        result_HS_charge(:,IES_no) = x(time * 6 + 1 + (IES_no - 1)* var :time * 7 +  (IES_no - 1)* var);
        eval(['EH',num2str(IES_no), '_Edr =  x(time * 7 + 1 + (IES_no - 1)* var :time * 8 +  (IES_no - 1)* var);']);
        eval(['EH',num2str(IES_no), '_Hdr =  x(time * 8 + 1 + (IES_no - 1)* var :time * 9 +  (IES_no - 1)* var);']);
    end
    for IES_no = 1 : number
        eval(['ES_para = ES',num2str(IES_no),'_para;']);
        eval(['HS_para = HS',num2str(IES_no),'_para;']);
        for pt = 1 : time
            result_ES_SOC(pt + 1, IES_no) = result_ES_SOC(pt, IES_no) - result_ES_discharge(pt, IES_no) / ES_para(7) / ES_para(1)...
                + result_ES_charge(pt, IES_no) * ES_para(7) / ES_para(1);
            result_HS_SOC(pt + 1, IES_no) = result_HS_SOC(pt, IES_no) - result_HS_discharge(pt, IES_no) / HS_para(7) / HS_para(1)...
                + result_HS_charge(pt, IES_no) * HS_para(7) / HS_para(1);
        end
    end
elseif isCentral == 2
    [result_Ele(:,1), result_CHP_G(:,1), result_Boiler_G(:,1), result_ES_discharge(:,1), result_ES_charge(:,1), result_HS_discharge(:,1), result_HS_charge(:,1), result_ES_SOC(:,1), result_HS_SOC(:,1), EH1_Le, EH1_Lh, EH1_solarP, EH1_windP, EH1_Edr, EH1_Hdr] = EH1.getResult;
    [result_Ele(:,2), result_CHP_G(:,2), result_Boiler_G(:,2), result_ES_discharge(:,2), result_ES_charge(:,2), result_HS_discharge(:,2), result_HS_charge(:,2), result_ES_SOC(:,2), result_HS_SOC(:,2), EH2_Le, EH2_Lh, EH2_solarP, EH2_windP, EH2_Edr, EH2_Hdr] = EH2.getResult;
    [result_Ele(:,3), result_CHP_G(:,3), result_Boiler_G(:,3), result_ES_discharge(:,3), result_ES_charge(:,3), result_HS_discharge(:,3), result_HS_charge(:,3), result_ES_SOC(:,3), result_HS_SOC(:,3), EH3_Le, EH3_Lh, EH3_solarP, EH3_windP, EH3_Edr, EH3_Hdr] = EH3.getResult;
end
main_handle_171013_v2
