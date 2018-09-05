% 集中式优化：单次
     para_init;

    [EH1_f, EH1_ub, EH1_lb, EH1_Aeq, EH1_beq, EH1_A, EH1_b, EH1_A_eleLimit_total ] = OptMatrix(eleLimit1, gasLimit1, EH1_Le, EH1_Lh,EH1_solarP, EH1_windP, CHP1_para, Boiler1_para, ES1_para, HS1_para,...
        gasPrice1, EH1_Le_drP_rate, EH1_Le_drP_total, EH1_Lh_drP_rate, EH1_Lh_drP_total);
    [EH2_f, EH2_ub, EH2_lb, EH2_Aeq, EH2_beq, EH2_A, EH2_b, EH2_A_eleLimit_total ] = OptMatrix(eleLimit2, gasLimit2, EH2_Le, EH2_Lh,EH2_solarP, EH2_windP, CHP2_para, Boiler2_para, ES2_para, HS2_para,...
        gasPrice1, EH2_Le_drP_rate, EH2_Le_drP_total, EH2_Lh_drP_rate, EH2_Lh_drP_total);
    [EH3_f, EH3_ub, EH3_lb, EH3_Aeq, EH3_beq, EH3_A, EH3_b, EH3_A_eleLimit_total ] = OptMatrix(eleLimit3, gasLimit3, EH3_Le, EH3_Lh,EH3_solarP, EH3_windP, CHP3_para, Boiler3_para, ES3_para, HS3_para,...
        gasPrice3, EH3_Le_drP_rate, EH3_Le_drP_total, EH3_Lh_drP_rate, EH3_Lh_drP_total);
    time = 24; %总时间段
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
    isCentral = 1;
    main_handle_171013_v2