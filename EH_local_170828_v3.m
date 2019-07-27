% clc
% clear

%多能协调优化
%能源枢纽的本地集中式优化策略
%2017.8.22 ver1.0
%2017.8.28 ver1.1 电、热功率平衡改为大于约束，增加购气量和的约束，BUG为解决：热有盈余的话储热会又充又放
%2017.8.28 ver2.0 改为class
%2018.1.10 ver3.0 market信息改用全局变量，对储能的又充又放采用等效的线性优化方法，时间间隔改为15min


classdef EH_local_170828_v3 < handle
    properties %可以设初始值
        %负荷与可再生能源
        Le_real;
        Le_pre;
        Lh_real;
        Lh_pre;
        solarP_real;
        solarP_pre;
        solarP_rate;
        windP_real;
        windP_rate;
        windP_pre;
        dev_L; %百分数
        dev_PV;
        dev_WT;
        
        %新增dr
        Le_drP_rate;
        Le_drP_total;
        Lh_drP_rate;
        Lh_drP_total;
        Le_T;
        Lh_T;
        
        %CHP的参数
        CHP_GE_eff;
        CHP_GH_eff;
        CHP_G_max;
        CHP_G_min;
        CHP_G_ramp
        
        %锅炉
        Boiler_eff;
        Boiler_G_max;
        % Boiler_G_min = 0;
        
        %联络线电交换和气交换约束
        Ele_max;
        Ele_min;
        Ele_eff;
        Gas_max;
        % Gas_min = 0;
        
        %电储能
        ES_totalC;
        ES_maxSOC;
        ES_minSOC;
        % ES_currentC; %随时间变化的
        ES_targetSOC;
        ES_Pmax;
        ES_eff;
        ES_selfd; % 1-自放电率
        %热储能
        HS_totalC;
        HS_maxSOC;
        HS_minSOC;
        % HS_currentC; %随时间变化的
        HS_targetSOC;
        HS_Hmax;
        HS_eff;
        HS_selfd; % 1-自放热率
        %投标
        demand_curve;
        
        %最终优化结果
        ES_SOC;
        HS_SOC;
        result_Ele;
        result_CHP_G;
        result_Boiler_G;
        result_ES_discharge;
        result_ES_charge;
        result_HS_discharge;
        result_HS_charge;
        result_E_dr;
        result_H_dr;
    end
    
    methods
        
        function obj = EH_local_170828_v3(eleLimit, gasLimit, Le_base, Lh_base, solar_base, wind_base, ...
                CHP_para, Boiler_para, ES_para, HS_para, dev_load, dev_solar, dev_wind, solar_rate, wind_rate,...
                Le_dr_rate, Le_dr_total, Lh_dr_rate, Lh_dr_total, Le_T, Lh_T)
            
            global priceNumbers period
            
            %联络线电交换和气交换约束
            obj.Ele_max = eleLimit(1);
            obj.Ele_min = eleLimit(2);
            obj.Ele_eff = eleLimit(3); %购电效率，即线路传输的线损率
            obj.Gas_max = gasLimit; %两个购气量的和
            % Gas_min = 0;
            
            %负荷与可再生能源
            obj.Le_real = Le_base; %base是日前初步的预测曲线
            obj.Lh_real = Lh_base;
            obj.solarP_real = solar_base;
            obj.solarP_rate = solar_rate;
            obj.windP_real = wind_base;
            obj.windP_rate = wind_rate;
            obj.dev_L = dev_load; %百分数
            obj.dev_PV = dev_solar;
            obj.dev_WT = dev_wind;
            
            %新增dr参数
            obj.Le_drP_rate = Le_dr_rate;
            obj.Le_drP_total = Le_dr_total;
            obj.Lh_drP_rate = Lh_dr_rate;
            obj.Lh_drP_total = Lh_dr_total;
            obj.Le_T = Le_T * obj.Le_drP_rate;
            obj.Lh_T = Lh_T * obj.Lh_drP_rate;
            
            %CHP的参数
            obj.CHP_GE_eff = CHP_para(1);
            obj.CHP_GH_eff = CHP_para(2);
            obj.CHP_G_max = CHP_para(3) / obj.CHP_GE_eff; %由额定电功率得到最大购气量
            obj.CHP_G_min = CHP_para(4) * obj.CHP_G_max;
            obj.CHP_G_ramp = CHP_para(5) * obj.CHP_G_max;
            
            %锅炉
            obj.Boiler_eff = Boiler_para(1);
            obj.Boiler_G_max = Boiler_para(2) / obj.Boiler_eff;  %10
            
            %电储能和热储能
            obj.ES_totalC = ES_para(1);
            obj.ES_maxSOC = ES_para(2);
            obj.ES_minSOC = ES_para(3);
            obj.ES_targetSOC = ES_para(4);
            obj.ES_Pmax = obj.ES_totalC / ES_para(6);
            obj.ES_eff = ES_para(7);
            obj.ES_selfd = 0.95;
            
            obj.HS_totalC = HS_para(1);
            obj.HS_maxSOC = HS_para(2);
            obj.HS_minSOC = HS_para(3);
            obj.HS_targetSOC = HS_para(4);
            obj.HS_Hmax = obj.HS_totalC / HS_para(6);
            obj.HS_eff = HS_para(7);
            obj.HS_selfd = 0.9;
            % 投标
            obj.demand_curve = zeros(priceNumbers, 1); %初始化投标曲线
            
            %最终优化结果
            obj.ES_SOC = zeros(24*period+1,1);
            obj.ES_SOC(1) = ES_para(4);
            obj.HS_SOC = zeros(24*period+1,1);
            obj.HS_SOC(1) = HS_para(4);
            obj.result_Ele = zeros(24*period,1);
            
            %新增可平移负荷水平
            obj.result_E_dr = zeros (24*period,1);
            obj.result_H_dr = zeros (24*period,1);
             
            obj.result_CHP_G = zeros(24*period,1);
            obj.result_Boiler_G = zeros(24*period,1);
            obj.result_ES_discharge = zeros(24*period,1);
            obj.result_ES_charge = zeros(24*period,1);
            obj.result_HS_discharge = zeros(24*period,1);
            obj.result_HS_charge = zeros(24*period,1);
        end
        
        % 可再生能源与负荷的预测业务，一天预测25次，日前一次，日内每小时一次
        function predict(obj, t_current) % 第几小时的预测，time=0是日前，1-24是日内
            global period
            if t_current ~= 0
                %rand 生成均匀分布的伪随机数 分布在（0~1）之间
                %randn 生成标准正态分布的伪随机数 （均值为0，方差为1）
                %                 randn('seed', t_current);
                
                %电、热负荷
                %                 Le_error = zeros(24*period,1);
                %                 Lh_error = zeros(24*period,1);
                %                 Le_error(t_current:24*period) = randn([(24*period-t_current+1),1]) .* obj.Le(t_current:24*period) * obj.dev_L; %预测误差
                %                 Lh_error(t_current:24*period) = randn([(24*period-t_current+1),1]) .* obj.Lh(t_current:24*period) * obj.dev_L;
                %                 obj.Le = obj.Le + Le_error; %更新
                %                 obj.Lh = obj.Lh + Lh_error;
                %改为只计当前时刻的预测误差，否则累计的误差过大；频繁重复预测也不科学
                for pt = t_current: 24 * period
                    if pt == t_current
                        dev_l = 0;
                        dev_res = 0;
                    else
                        dev_l = 0.08;
                        dev_res = 0.1;
                    end
                    Le_error = randn() * obj.Le_real(pt) * dev_l; %预测误差
                    Lh_error = randn() * obj.Lh_real(pt) * dev_l;
                    obj.Le_pre(pt) = obj.Le_real(pt) + Le_error; %更新
                    obj.Lh_pre(pt) = obj.Lh_real(pt) + Lh_error;
                    solarP_error = randn() * obj.solarP_real(pt) * dev_res; %预测误差，solarP(t_current)=0的时候自然为零
                    windP_error = randn() * obj.windP_real(pt) * dev_res; %预测误差             
                    obj.solarP_pre(pt) = obj.solarP_real(pt) + solarP_error; %更新  
                    obj.windP_pre(pt) = obj.windP_real(pt) + windP_error;
                end                
            else
                for pt = 1: 24 * period
                    dev_l = 0.2;
                    dev_res = 0.3;                    
                    Le_error = randn() * obj.Le_real(pt) * dev_l; %预测误差
                    Lh_error = randn() * obj.Lh_real(pt) * dev_l;
                    obj.Le_pre(pt) = obj.Le_real(pt) + Le_error; %更新
                    obj.Lh_pre(pt) = obj.Lh_real(pt) + Lh_error;
                    solarP_error = randn() * obj.solarP_real(pt) * dev_res; %预测误差，solarP(t_current)=0的时候自然为零
                    windP_error = randn() * obj.windP_real(pt) * dev_res; %预测误差             
                    obj.solarP_pre(pt) = obj.solarP_real(pt) + solarP_error; %更新  
                    obj.windP_pre(pt) = obj.windP_real(pt) + windP_error;
                end                    
            end
            obj.Le_pre(obj.Le_pre < 0) = 0;
            obj.Lh_pre(obj.Lh_pre < 0) = 0;
            obj.solarP_pre(obj.solarP_pre < 0) = 0;
            obj.solarP_pre(obj.solarP_pre > obj.solarP_rate) = obj.solarP_rate;
            obj.windP_pre(obj.windP_pre < 0) = 0;
            obj.windP_pre(obj.windP_pre > obj.windP_rate) = obj.windP_rate;
        end
        
        %接收市场出清价格，做自治优化，更新自身状态
        function [x,fval,exitflag,output,lambda] = handlePrice(obj, Eprice, Gprice, t_current) %这里的x随着t_current会越来越少
            [x,fval,exitflag,output,lambda] = localOptimal(obj, Eprice, Gprice, t_current, 9e9); % conditionEle = 9e9
        end
        
        
        %用于TC单次出清
        %接收市场出清价格，按投标曲线得到出清功率，在保持出清功率不变的情况下，求解最优化
        %{
        function [x,fval,exitflag,output,lambda] = conditionHandlePrice(obj, Eprice, Gprice, t_current) %这里的x随着t_current会越来越少
            global period
            conditionEle = getClearDemand(obj.demand_curve, Eprice(t_current));
            
            [x,fval,exitflag,output,lambda] = localOptimal(obj, Eprice, Gprice, t_current, conditionEle);
            time = 24*period - t_current + 1; %总时间段
            
            % 只执行当前周期的结果
            obj.result_Ele(t_current) = x(1);
            obj.result_CHP_G(t_current) = x(time+1);
            obj.result_Boiler_G(t_current) = x(time*2+1);
            obj.result_ES_discharge(t_current) = x(time*3+1);
            obj.result_ES_charge(t_current) = x(time*4+1);
            obj.result_HS_discharge(t_current) = x(time*5+1);
            obj.result_HS_charge(t_current) = x(time*6+1);
            
            %更新储能状态
            obj.ES_SOC(t_current+1) = obj.ES_SOC(t_current) - obj.result_ES_discharge(t_current) / obj.ES_eff / obj.ES_totalC + obj.result_ES_charge(t_current) * obj.ES_eff / obj.ES_totalC;
            obj.HS_SOC(t_current+1) = obj.HS_SOC(t_current) - obj.result_HS_discharge(t_current) / obj.HS_eff / obj.HS_totalC + obj.result_HS_charge(t_current) * obj.HS_eff / obj.HS_totalC;
        end
        %}
        
        % 接收市场出清价格、出清功率，在保持出清功率不变的情况下，求解最优化
        % 只更新当前状态
        function [x,fval,exitflag,output,lambda] = conditionHandlePrice_2(obj, Eprice, Gprice, t_current, clearDemand) % 这里的x随着t_current会越来越少
            global period
            conditionEle = clearDemand;
            
            [x,fval,exitflag,output,lambda] = localOptimal(obj, Eprice, Gprice, t_current, conditionEle);
            time = 24*period - t_current + 1; %总时间段
            
            % 只执行当前周期的结果
            obj.result_Ele(t_current) = x(1);
            obj.result_CHP_G(t_current) = x(time+1);
            obj.result_Boiler_G(t_current) = x(time*2+1);
            obj.result_ES_discharge(t_current) = x(time*3+1);
            obj.result_ES_charge(t_current) = x(time*4+1);
            obj.result_HS_discharge(t_current) = x(time*5+1);
            obj.result_HS_charge(t_current) = x(time*6+1);
            obj.result_E_dr(t_current) = x(time*7+1);
            obj.result_H_dr(t_current) = x(time*8+1);
            %更新储能状态
            obj.ES_SOC(t_current+1) = obj.ES_selfd * obj.ES_SOC(t_current) - obj.result_ES_discharge(t_current) / obj.ES_eff / obj.ES_totalC + obj.result_ES_charge(t_current) * obj.ES_eff / obj.ES_totalC;
            obj.HS_SOC(t_current+1) = obj.HS_selfd * obj.HS_SOC(t_current) - obj.result_HS_discharge(t_current) / obj.HS_eff / obj.HS_totalC + obj.result_HS_charge(t_current) * obj.HS_eff / obj.HS_totalC;
        end
        
        
        % _DA为全时段日前优化，因此一次性更新后续所有时段的信息
        function [x,fval,exitflag,output,lambda] = conditionHandlePrice_all(obj, Eprice, Gprice, t_current, clearDemand) %这里的x随着t_current会越来越少
            global period
            conditionEle = clearDemand;
            
            [x,fval,exitflag,output,lambda] = localOptimal(obj, Eprice, Gprice, t_current, conditionEle);
            time = 24*period - t_current + 1; %总时间段
            for pt = t_current: 24 * period
                obj.result_Ele(pt) = x(1 + pt - t_current);
                obj.result_CHP_G(pt) = x(time + 1 + pt - t_current);
                obj.result_Boiler_G(pt) = x(time*2 + 1 + pt - t_current);
                obj.result_ES_discharge(pt) = x(time*3 + 1 + pt - t_current);
                obj.result_ES_charge(pt) = x(time*4 + 1 + pt - t_current);
                obj.result_HS_discharge(pt) = x(time*5 + 1 + pt - t_current);
                obj.result_HS_charge(pt) = x(time*6 + 1 + pt - t_current);
                %新增可平移负荷水平
                obj.result_E_dr(pt) = x(time*7+ 1 +pt -t_current);
                obj.result_H_dr(pt) = x(time*8+ 1 +pt -t_current);
                %更新储能状态
                obj.ES_SOC(pt+1) = obj.ES_selfd * obj.ES_SOC(pt) - obj.result_ES_discharge(pt) / obj.ES_eff / obj.ES_totalC + obj.result_ES_charge(pt) * obj.ES_eff / obj.ES_totalC;
                obj.HS_SOC(pt+1) = obj.HS_selfd * obj.HS_SOC(pt) - obj.result_HS_discharge(pt) / obj.HS_eff / obj.HS_totalC + obj.result_HS_charge(pt) * obj.HS_eff / obj.HS_totalC;
                
            end
        end
        % 输出优化结果
        function [Ele, G_CHP, G_Boiler, ES_discharge, ES_charge, HS_discharge, HS_charge, ES_SOC, HS_SOC, Le, Lh, solarP, windP, Edr, Hdr] = getResult(obj)
            Ele = obj.result_Ele;
            G_CHP = obj.result_CHP_G;
            G_Boiler = obj.result_Boiler_G;
            ES_discharge = obj.result_ES_discharge;
            ES_charge = obj.result_ES_charge;
            HS_discharge = obj.result_HS_discharge;
            HS_charge = obj.result_HS_charge;
            ES_SOC = obj.ES_SOC;
            HS_SOC = obj.HS_SOC;
            % 同时输出偏移后的负荷与可再生能源
            Le = obj.Le_real;
            Lh = obj.Lh_real;
            solarP = obj.solarP_real;
            windP = obj.windP_real;
            %新增负荷水平
            Edr= obj.result_E_dr;
            Hdr= obj.result_H_dr;
        end
        
        
        
        % 测试优化结果
        function [result_balance_P, result_balance_H, result_check_ES, result_check_HS] = testResult(obj)
            %电、热功率平衡性测试，降低要求至大于零
            result_balance_P = obj.result_Ele + obj.CHP_GE_eff.*obj.result_CHP_G + obj.result_ES_discharge - obj.result_ES_charge - obj.Le_pre + obj.windP_pre + obj.solarP_pre -obj.result_E_dr;
            result_balance_H = obj.CHP_GH_eff.*obj.result_CHP_G + obj.Boiler_eff.*obj.result_Boiler_G + obj.result_HS_discharge - obj.result_HS_charge - obj.Lh_pre-obj.result_H_dr;
            %充、放功率至少有一个是零
            result_check_ES = obj.result_ES_discharge .* obj.result_ES_charge;
            result_check_HS = obj.result_HS_discharge .* obj.result_HS_charge; %有一点小问题，因为热过于充裕
        end
        
        
        % 解最优化
        function [x,fval,exitflag,output,lambda] = localOptimal(obj, Eprice, Gprice, t_current, conditionEle)
            %参数里面必须要加obj，否则会报错too many input arguments
            
            %当conditionEle<9e9时，表示在给定的当前时刻的购电量下，解最优化
            %将当前时刻的购电量，最小值和最大值都设定为给定值
            
            global period
            time = 24 * period - t_current + 1; %总时间段
            var = time * 9; %总变量数
            %第1,2,3组time是购电量、CHP购气量、锅炉购气量，第4-7组time是储电、储热的放、充功率
            %第8-9组是当前电和热的可平移负荷水平
            
            %求一次项系数f 是个列向量
            f = zeros(var, 1);
            for i = 1 : time
                % f(i, 1) = Eprice(i); %Eprice的size是随t_current变化的
                f(i, 1) = Eprice(t_current + i - 1); %Eprice的size不变，但是取部分值
                f(time+i, 1) = Gprice;
                f(time*2+i, 1) = Gprice;
            end
            
            %变量上下限
            ub = zeros(var, 1);
            lb = zeros(var, 1);
            for i = 1 : time
                ub(i, 1) = obj.Ele_max;
                ub(time + i, 1) = obj.CHP_G_max;
                ub(time * 2 + i, 1) = obj.Boiler_G_max;
                ub(time * 3 + i, 1) = obj.ES_Pmax;
                ub(time * 4 + i, 1) = obj.ES_Pmax;
                ub(time * 5 + i, 1) = obj.HS_Hmax;    
                ub(time * 6 + i, 1) = obj.HS_Hmax;
            end
            ub(time * 7 + 1:time * 8, 1) = obj.Le_T(t_current: end);
            ub(time * 8 + 1:time * 9, 1) = obj.Lh_T(t_current: end);

            for i = 1 : time
                lb(i, 1) = obj.Ele_min;
                lb(time + i, 1) = obj.CHP_G_min;
            end
            if length(conditionEle)>1
                ub(1:length(conditionEle),1) = conditionEle;
                lb(1:length(conditionEle),1) = conditionEle;
            else
                if conditionEle < 9e9
                    ub(1,1) = conditionEle;
                    lb(1,1) = conditionEle;
                end
            end
            
            %等式约束包括：电、热平衡约束（供大于求，改为不等式）；电、热储能平衡性约束（改为不等式约束？）
            %电、热平衡约束
            Aeq_Ebus = zeros(time, var);
            Aeq_Hbus = zeros(time, var);
            beq_Ebus = - obj.Le_pre(t_current : 24*period) + obj.windP_pre(t_current : 24*period) + obj.solarP_pre(t_current : 24*period); %Le的size不变，但是取部分值
            beq_Hbus = - obj.Lh_pre(t_current : 24*period); %Lh的size不变，但是取部分值
            for i=1:time
                Aeq_Ebus(i, i) = - obj.Ele_eff; %线损率
                Aeq_Ebus(i, time + i) = - obj.CHP_GE_eff;
                Aeq_Ebus(i, time * 3 + i) = - 1; %放电
                Aeq_Ebus(i, time * 4 + i) = 1; %充电
                Aeq_Ebus(i, time * 7 + i) = 1;
            end
            for i=1:time
                Aeq_Hbus(i,time + i) = - obj.CHP_GH_eff;
                Aeq_Hbus(i,time * 2 + i) = - obj.Boiler_eff;
                Aeq_Hbus(i,time * 5 + i) = - 1; %放热
                Aeq_Hbus(i,time * 6 + i) = 1; %充热
                Aeq_Hbus(i,time * 8 + i) = 1;
            end
            
            %可平移负荷约束
            Aeq_Edr = zeros(1, var);
            beq_Edr = obj.Le_drP_total - sum(obj.result_E_dr(1:t_current-1));
            for i=1:time
                Aeq_Edr(1, time*7+i) = 1;
            end
            
            Aeq_Hdr = zeros(1, var);
            beq_Hdr = obj.Lh_drP_total - sum(obj.result_H_dr(1:t_current-1));
            for i=1:time
                Aeq_Hdr(1, time*8+i) = 1;
            end
            %电、热储能平衡性约束
            Aeq_ES = zeros(1, var);
            beq_ES = - (obj.ES_targetSOC - obj.ES_SOC(t_current) * obj.ES_selfd ^ time) * obj.ES_totalC;
            for i=1:time
                Aeq_ES(1, time * 3 + i) = obj.ES_selfd ^ (time - i) / obj.ES_eff; %放电
                Aeq_ES(1, time * 4 + i) = - obj.ES_selfd ^ (time - i) * obj.ES_eff; %充电
            end
            Aeq_HS = zeros(1, var);
            beq_HS = - (obj.HS_targetSOC - obj.HS_SOC(t_current) * obj.HS_selfd ^ time) * obj.HS_totalC;
            for i=1:time
                Aeq_HS(1, time * 5 + i) = obj.HS_selfd ^ (time - i) / obj.HS_eff; %放热
                Aeq_HS(1, time * 6 + i) = - obj.HS_selfd ^ (time - i) * obj.HS_eff; %充热
            end
            %不等式约束包括：SOC约束；购气量和的约束；（爬坡率约束）
            %SOC约束 A1是上限，A2是下限
            A1_Esoc = zeros(time, var);
            b1_Esoc = zeros(time, 1);
            b2_Esoc = zeros(time, 1);
            for i = 1 : time
               b1_Esoc(i, 1) = (obj.ES_maxSOC - obj.ES_SOC(t_current) * obj.ES_selfd ^ i) * obj.ES_totalC; 
               b2_Esoc(i, 1) = (obj.ES_selfd ^ i * obj.ES_SOC(t_current) - obj.ES_minSOC) * obj.ES_totalC;
            end
            
            for i=1:time
                for j=1 : i
                    A1_Esoc(i, time * 3 + j) = -  obj.ES_selfd ^ (i - j) / obj.ES_eff; %放电
                    A1_Esoc(i, time * 4 + j) = obj.ES_selfd ^ (i - j) * obj.ES_eff; %充电
                end
            end
            A2_Esoc = - A1_Esoc;
            A1_Hsoc = zeros(time, var);
            b1_Hsoc = zeros(time, 1);
            b2_Hsoc = zeros(time, 1);
            for i = 1 : time
                b1_Hsoc(i, 1) = (obj.HS_maxSOC - obj.HS_SOC(t_current) * obj.HS_selfd ^ i) * obj.HS_totalC;
                b2_Hsoc(i, 1) = (obj.HS_SOC(t_current) * obj.HS_selfd ^ i - obj.HS_minSOC) * obj.HS_totalC;
            end
            for i=1:time
                for j=1 : i
                    A1_Hsoc(i, time*5+j) = - obj.HS_selfd ^ (i - j) / obj.HS_eff; %放热
                    A1_Hsoc(i, time*6+j) = obj.HS_selfd ^ (i - j) * obj.HS_eff; %充热
                end
            end
            A2_Hsoc = -A1_Hsoc;
            %购气量和的约束
            A_Gmax = zeros(time, var);
            b_Gmax = ones(time,1) .* obj.Gas_max;
            for i=1:time
                A_Gmax(i, time+i) = 1;
                A_Gmax(i, time*2+i) = 1;
            end
            
            %CHP机组爬坡率约束
            A_ramp_tmp = diag(ones(1, time))-diag(ones(1, time -1),1);
            A_ramp_tmp = A_ramp_tmp(1: end -1,:);
            A_ramp = zeros(time -1, var); 
            A_ramp(time * 1 + 1 : time * 2, :) = A_ramp_tmp;
            b_ramp = obj.CHP_G_ramp * ones(time-1, 1);
            %归纳所有线性约束
            %等式约束包括：电、热平衡约束（改为不等式），电、热储能平衡性约束（改为不等式）
            %不等式约束包括：SOC约束，购气量和的约束，（爬坡率约束）
            Aeq=[Aeq_Edr; Aeq_Hdr; Aeq_Hbus; Aeq_ES; Aeq_HS;];
            beq=[beq_Edr; beq_Hdr; beq_Hbus';beq_ES; beq_HS;];
            A=[ Aeq_Ebus; A1_Esoc; A2_Esoc; A1_Hsoc; A2_Hsoc; A_Gmax; A_ramp; -A_ramp];
            b=[ beq_Ebus'; b1_Esoc; b2_Esoc; b1_Hsoc; b2_Hsoc; b_Gmax; b_ramp; b_ramp];
 
            [x,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub);
            
            %             options = optimoptions('fmincon','MaxFunEvals',1000000);
            %             [x,fval,exitflag,output,lambda] = fmincon('myfun_1', x0, A, b, Aeq, beq, lb, ub, 'mycon', options);
            
            if exitflag ~= 1
                %                 error('没有可行解')
            end
            
        end
        
        function [f, intcon, A, b, Aeq, beq, lb, ub, A_eleLimit_total] = MilpMatrix(obj, Gprice, t_current)
            global period elePrice
            time = 24 * period - t_current + 1;
            var = time * 13; %总变量数
            intcon = [time * 9 + 1 : time * 13];%charge/dischage binary
            
            %求一次项系数f 是个列向量
            f = zeros(var, 1);
            for i = 1 : time
                f(i, 1) = elePrice(t_current + i - 1); %Eprice的size不变，但是取部分值
                f(time + i, 1) = Gprice;
                f(time * 2 + i, 1) = Gprice;
            end
            
            %变量上下限
            ub = zeros(var, 1);
            lb = zeros(var, 1);
            for i = 1 : time
                ub(i, 1) = obj.Ele_max;
                ub(time + i, 1) = obj.CHP_G_max;
                ub(time * 2 + i, 1) = obj.Boiler_G_max;
                ub(time * 3 + i, 1) = obj.ES_Pmax;
                ub(time * 4 + i, 1) = obj.ES_Pmax;
                ub(time * 5 + i, 1) = obj.HS_Hmax;
                ub(time * 6 + i, 1) = obj.HS_Hmax;
                ub(time * 9 + i, 1) = 1;
                ub(time * 10 + i, 1) = 1;
                ub(time * 11 + i, 1) = 1;
                ub(time * 12 + i, 1) = 1;
            end
            ub(time * 7 + 1:time * 8, 1) = obj.Le_T(t_current: end);
            ub(time * 8 + 1:time * 9, 1) = obj.Lh_T(t_current: end);
            
            for i = 1 : time
                lb(i, 1) = obj.Ele_min;
                lb(time + i, 1) = obj.CHP_G_min;
            end
            
            %等式约束包括：电、热平衡约束（供大于求，改为不等式）；电、热储能平衡性约束（改为不等式约束？）
            %电、热平衡约束
            Aeq_Ebus = zeros(time, var);
            Aeq_Hbus = zeros(time, var);
            beq_Ebus = - obj.Le_pre(t_current : 24 * period) + obj.windP_pre(t_current : 24 * period) + obj.solarP_pre(t_current : 24 * period);
            beq_Hbus = - obj.Lh_pre(t_current : 24 * period); 
            for i=1:time
                Aeq_Ebus(i,i) = - obj.Ele_eff; %线损率
                Aeq_Ebus(i,time + i) = - obj.CHP_GE_eff;
                Aeq_Ebus(i,time * 3 + i) = - 1; %放电
                Aeq_Ebus(i,time * 4 + i) = 1; %充电
                Aeq_Ebus(i,time * 7 + i) = 1;
            end
            for i=1:time
                Aeq_Hbus(i,time + i) = - obj.CHP_GH_eff;
                Aeq_Hbus(i,time * 2 + i) = - obj.Boiler_eff;
                Aeq_Hbus(i,time * 5 + i) = - 1; %放热
                Aeq_Hbus(i,time * 6 + i) = 1; %充热
                Aeq_Hbus(i,time * 8 + i) = 1;
            end
            
            %可平移负荷约束
            Aeq_Edr = zeros(1, var);
            beq_Edr = obj.Le_drP_total - sum(obj.result_E_dr(1:t_current-1));
            for i=1:time
                Aeq_Edr(1, time * 7 + i) = 1;
            end
            
            Aeq_Hdr = zeros(1, var);
            beq_Hdr = obj.Lh_drP_total - sum(obj.result_H_dr(1:t_current-1));
            for i=1:time
                Aeq_Hdr(1, time * 8 + i) = 1;
            end
            %电、热储能平衡性约束
            Aeq_ES = zeros(1, var);
            beq_ES = - (obj.ES_targetSOC - obj.ES_SOC(t_current) * obj.ES_selfd ^ time) * obj.ES_totalC;
            for i=1:time
                Aeq_ES(1, time * 3 + i) = obj.ES_selfd ^ (time - i) / obj.ES_eff; %放电
                Aeq_ES(1, time * 4 + i) = - obj.ES_selfd ^ (time - i) * obj.ES_eff; %充电
            end
            Aeq_HS = zeros(1, var);
            beq_HS = - (obj.HS_targetSOC - obj.HS_SOC(t_current) * obj.HS_selfd ^ time) * obj.HS_totalC;
            for i=1:time
                Aeq_HS(1, time * 5 + i) = obj.HS_selfd ^ (time - i) / obj.HS_eff; %放热
                Aeq_HS(1, time * 6 + i) = - obj.HS_selfd ^ (time - i) * obj.HS_eff; %充热
            end
            
            %不等式约束
            %EES binary约束
            A_EESbinary = zeros(time, var);
            for i=1:time
                A_EESbinary(i, time * 9 + i) = 1; %放电
                A_EESbinary(i, time * 10 + i) = 1; %充电
            end
            b_EESbinary = ones(time, 1);
            
            %充放电上下限功率
            A_EEScharge = zeros(time, var);
            for i=1:time
                A_EEScharge(i, time * 10 + i) = -obj.ES_Pmax; 
                A_EEScharge(i,time * 4 + i) = 1; 

            end
            b_EEScharge = ones(time, 1);
            
            A_EESdischarge = zeros(time, var);
            for i=1:time
                A_EESdischarge(i, time * 9 + i) = -obj.ES_Pmax; 
                A_EESdischarge(i,time * 3 + i) = 1; 

            end
            b_EESdischarge = ones(time, 1);
            
            %TES binary约束
            A_TESbinary = zeros(time, var);
            for i=1:time
                A_TESbinary(i, time * 11 + i) = 1; %放电
                A_TESbinary(i, time * 12 + i) = 1; %充电
            end
            b_TESbinary = ones(time, 1);
            
            %充放电上下限功率
            A_TEScharge = zeros(time, var);
            for i=1:time
                A_TEScharge(i, time * 12 + i) = -obj.HS_Hmax; 
                A_TEScharge(i,time * 6 + i) = 1; 

            end
            b_TEScharge = ones(time, 1);
            
            A_TESdischarge = zeros(time, var);
            for i=1:time
                A_TESdischarge(i, time * 11 + i) = -obj.HS_Hmax; 
                A_TESdischarge(i,time * 5 + i) = 1; 

            end
            b_TESdischarge = ones(time, 1);
            
            A_binary = [A_EESbinary; A_TESbinary; A_EEScharge; A_EESdischarge; A_TEScharge; A_TESdischarge];
            b_binary = [b_EESbinary; b_TESbinary; b_EEScharge; b_EESdischarge; b_TEScharge; b_TESdischarge];
            
              %不等式约束包括：SOC约束；购气量和的约束；（爬坡率约束）
            %SOC约束 A1是上限，A2是下限
            A1_Esoc = zeros(time, var);
            b1_Esoc = zeros(time, 1);
            b2_Esoc = zeros(time, 1);
            for i = 1 : time
               b1_Esoc(i, 1) = (obj.ES_maxSOC - obj.ES_SOC(t_current) * obj.ES_selfd ^ i) * obj.ES_totalC; 
               b2_Esoc(i, 1) = (obj.ES_selfd ^ i * obj.ES_SOC(t_current) - obj.ES_minSOC) * obj.ES_totalC;
            end
            
            for i=1:time
                for j=1 : i
                    A1_Esoc(i, time * 3 + j) = -  obj.ES_selfd ^ (i - j) / obj.ES_eff; %放电
                    A1_Esoc(i, time * 4 + j) = obj.ES_selfd ^ (i - j) * obj.ES_eff; %充电
                end
            end
            A2_Esoc = - A1_Esoc;
            A1_Hsoc = zeros(time, var);
            b1_Hsoc = zeros(time, 1);
            b2_Hsoc = zeros(time, 1);
            for i = 1 : time
                b1_Hsoc(i, 1) = (obj.HS_maxSOC - obj.HS_SOC(t_current) * obj.HS_selfd ^ i) * obj.HS_totalC;
                b2_Hsoc(i, 1) = (obj.HS_SOC(t_current) * obj.HS_selfd ^ i - obj.HS_minSOC) * obj.HS_totalC;
            end
            for i=1:time
                for j=1 : i
                    A1_Hsoc(i, time*5+j) = - obj.HS_selfd ^ (i - j) / obj.HS_eff; %放热
                    A1_Hsoc(i, time*6+j) = obj.HS_selfd ^ (i - j) * obj.HS_eff; %充热
                end
            end
            A2_Hsoc = -A1_Hsoc;
            %购气量和的约束
            A_Gmax = zeros(time, var);
            b_Gmax = ones(time,1) .* obj.Gas_max;
            for i=1:time
                A_Gmax(i, time + i) = 1;
                A_Gmax(i, time * 2 + i) = 1;
            end
            %CHP机组爬坡率约束
            A_ramp_tmp = diag(ones(1, time))-diag(ones(1, time -1),1);
            A_ramp_tmp = A_ramp_tmp(1: end -1,:);
            A_ramp = zeros(time -1, var); 
            A_ramp(time * 1 + 1 : time * 2, :) = A_ramp_tmp;
            b_ramp = obj.CHP_G_ramp * ones(time-1, 1);
           
            %归纳所有线性约束
            %等式约束包括：电、热平衡约束（改为不等式），电、热储能平衡性约束（改为不等式）
            %不等式约束包括：SOC约束，购气量和的约束，（爬坡率约束）
            Aeq=[Aeq_Edr; Aeq_Hdr;Aeq_Hbus; Aeq_ES; Aeq_HS; ];
            beq=[beq_Edr; beq_Hdr;beq_Hbus';beq_ES; beq_HS;];
            A=[Aeq_Ebus; A1_Esoc; A2_Esoc; A1_Hsoc; A2_Hsoc; A_Gmax; A_binary; A_ramp; -A_ramp];
            b=[beq_Ebus'; b1_Esoc; b2_Esoc; b1_Hsoc; b2_Hsoc; b_Gmax; b_binary; b_ramp; b_ramp];

           % 需要额外增加一个购电量的上、下限约束
            A_eleLimit_total = zeros(time, var);
            for i=1:time
                A_eleLimit_total(i, i) = 1;
            end 
        end
        
        
        function [f, ub, lb, Aeq, beq, A, b, A_eleLimit_total] = OptMatrix(obj, Gprice, t_current)
            
            global period elePrice
            time = 24 * period - t_current + 1;
            var = time * 9; %总变量数

            %求一次项系数f 是个列向量
            f = zeros(var, 1);
            for i = 1 : time
                f(i, 1) = elePrice(t_current + i - 1); %Eprice的size不变，但是取部分值
                f(time + i, 1) = Gprice;
                f(time * 2 + i, 1) = Gprice;
            end
            
            %变量上下限
            ub = zeros(var, 1);
            lb = zeros(var, 1);
            for i = 1 : time
                ub(i, 1) = obj.Ele_max;
                ub(time + i, 1) = obj.CHP_G_max;
                ub(time * 2 + i, 1) = obj.Boiler_G_max;
                ub(time * 3 + i, 1) = obj.ES_Pmax;
                ub(time * 4 + i, 1) = obj.ES_Pmax;
                ub(time * 5 + i, 1) = obj.HS_Hmax;
                ub(time * 6 + i, 1) = obj.HS_Hmax; 
            end
            ub(time * 7 + 1:time * 8, 1) = obj.Le_T(t_current: end);
            ub(time * 8 + 1:time * 9, 1) = obj.Lh_T(t_current: end);
            for i = 1 : time
                lb(i, 1) = obj.Ele_min;
                lb(time + i, 1) = obj.CHP_G_min;
            end
            
            %等式约束包括：电、热平衡约束（供大于求，改为不等式）；电、热储能平衡性约束（改为不等式约束？）
            %电、热平衡约束
            Aeq_Ebus = zeros(time, var);
            Aeq_Hbus = zeros(time, var);
            beq_Ebus = - obj.Le_pre(t_current : 24 * period) + obj.windP_pre(t_current : 24 * period) + obj.solarP_pre(t_current : 24 * period);
            beq_Hbus = - obj.Lh_pre(t_current : 24 * period); 
            for i=1:time
                Aeq_Ebus(i,i) = - obj.Ele_eff; %线损率
                Aeq_Ebus(i,time + i) = - obj.CHP_GE_eff;
                Aeq_Ebus(i,time * 3 + i) = - 1; %放电
                Aeq_Ebus(i,time * 4 + i) = 1; %充电
                Aeq_Ebus(i,time * 7 + i) = 1;
            end
            for i=1:time
                Aeq_Hbus(i,time + i) = - obj.CHP_GH_eff;
                Aeq_Hbus(i,time * 2 + i) = - obj.Boiler_eff;
                Aeq_Hbus(i,time * 5 + i) = - 1; %放热
                Aeq_Hbus(i,time * 6 + i) = 1; %充热
                Aeq_Hbus(i,time * 8 + i) = 1;
            end
            
            %可平移负荷约束
            Aeq_Edr = zeros(1, var);
            beq_Edr = obj.Le_drP_total - sum(obj.result_E_dr(1:t_current-1));
            for i=1:time
                Aeq_Edr(1, time * 7 + i) = 1;
            end
            
            Aeq_Hdr = zeros(1, var);
            beq_Hdr = obj.Lh_drP_total - sum(obj.result_H_dr(1:t_current-1));
            for i=1:time
                Aeq_Hdr(1, time * 8 + i) = 1;
            end
            %电、热储能平衡性约束
            Aeq_ES = zeros(1, var);
            beq_ES = - (obj.ES_targetSOC - obj.ES_SOC(t_current) * obj.ES_selfd ^ time) * obj.ES_totalC;
            for i=1:time
                Aeq_ES(1, time * 3 + i) = obj.ES_selfd ^ (time - i) / obj.ES_eff; %放电
                Aeq_ES(1, time * 4 + i) = - obj.ES_selfd ^ (time - i) * obj.ES_eff; %充电
            end
            Aeq_HS = zeros(1, var);
            beq_HS = - (obj.HS_targetSOC - obj.HS_SOC(t_current) * obj.HS_selfd ^ time) * obj.HS_totalC;
            for i=1:time
                Aeq_HS(1, time * 5 + i) = obj.HS_selfd ^ (time - i) / obj.HS_eff; %放热
                Aeq_HS(1, time * 6 + i) = - obj.HS_selfd ^ (time - i) * obj.HS_eff; %充热
            end
            
             %不等式约束包括：SOC约束；购气量和的约束；（爬坡率约束）
            %SOC约束 A1是上限，A2是下限
            A1_Esoc = zeros(time, var);
            b1_Esoc = zeros(time, 1);
            b2_Esoc = zeros(time, 1);
            for i = 1 : time
               b1_Esoc(i, 1) = (obj.ES_maxSOC - obj.ES_SOC(t_current) * obj.ES_selfd ^ i) * obj.ES_totalC; 
               b2_Esoc(i, 1) = (obj.ES_selfd ^ i * obj.ES_SOC(t_current) - obj.ES_minSOC) * obj.ES_totalC;
            end
            
            for i=1:time
                for j=1 : i
                    A1_Esoc(i, time * 3 + j) = -  obj.ES_selfd ^ (i - j) / obj.ES_eff; %放电
                    A1_Esoc(i, time * 4 + j) = obj.ES_selfd ^ (i - j) * obj.ES_eff; %充电
                end
            end
            A2_Esoc = - A1_Esoc;
            A1_Hsoc = zeros(time, var);
            b1_Hsoc = zeros(time, 1);
            b2_Hsoc = zeros(time, 1);
            for i = 1 : time
                b1_Hsoc(i, 1) = (obj.HS_maxSOC - obj.HS_SOC(t_current) * obj.HS_selfd ^ i) * obj.HS_totalC;
                b2_Hsoc(i, 1) = (obj.HS_SOC(t_current) * obj.HS_selfd ^ i - obj.HS_minSOC) * obj.HS_totalC;
            end
            for i=1:time
                for j=1 : i
                    A1_Hsoc(i, time*5+j) = - obj.HS_selfd ^ (i - j) / obj.HS_eff; %放热
                    A1_Hsoc(i, time*6+j) = obj.HS_selfd ^ (i - j) * obj.HS_eff; %充热
                end
            end
            A2_Hsoc = -A1_Hsoc;
            %购气量和的约束
            A_Gmax = zeros(time, var);
            b_Gmax = ones(time,1) .* obj.Gas_max;
            for i=1:time
                A_Gmax(i, time + i) = 1;
                A_Gmax(i, time * 2 + i) = 1;
            end
           
            %CHP机组爬坡率约束
            A_ramp_tmp = diag(ones(1, time))-diag(ones(1, time -1),1);
            A_ramp_tmp = A_ramp_tmp(1: end -1,:);
            A_ramp = zeros(time -1, var); 
            A_ramp(:, time * 1 + 1 : time * 2) = A_ramp_tmp;
            b_ramp = obj.CHP_G_ramp * ones(time-1, 1);
            
            %归纳所有线性约束
            %等式约束包括：电、热平衡约束（改为不等式），电、热储能平衡性约束（改为不等式）
            %不等式约束包括：SOC约束，购气量和的约束，（爬坡率约束）
            Aeq=[Aeq_Edr; Aeq_Hdr;Aeq_Hbus; Aeq_ES; Aeq_HS; ];
            beq=[beq_Edr; beq_Hdr;beq_Hbus';beq_ES; beq_HS;];
            A=[Aeq_Ebus;  A1_Esoc; A2_Esoc; A1_Hsoc; A2_Hsoc; A_Gmax; A_ramp; -A_ramp];
            b=[beq_Ebus'; b1_Esoc; b2_Esoc; b1_Hsoc; b2_Hsoc; b_Gmax; b_ramp; b_ramp];

           % 需要额外增加一个购电量的上、下限约束
            A_eleLimit_total = zeros(time, var);
            for i=1:time
                A_eleLimit_total(i, i) = 1;
            end 
        end
        
         function  update_central(obj, x, t_current, IES_no, isMILP) %这里的x随着t_current会越来越少
            global period 
            time = 24 * period - t_current + 1; %总时间段
            if nargin < 5%LP
                var = time * 9;
            else%MILP
                var = time * 13;
            end
            obj.result_Ele(t_current) = x(1 + (IES_no - 1) * var);
            obj.result_CHP_G(t_current) = x(time + 1 + (IES_no - 1) * var);
            obj.result_Boiler_G(t_current) = x(time * 2 + 1 + (IES_no - 1) * var);
            obj.result_ES_discharge(t_current) = x(time * 3 + 1 + (IES_no - 1) * var);
            obj.result_ES_charge(t_current) = x(time * 4 + 1 + (IES_no - 1) * var );
            obj.result_HS_discharge(t_current) = x(time * 5 + 1 + (IES_no - 1) * var);
            obj.result_HS_charge(t_current) = x(time * 6 + 1 + (IES_no - 1) * var);
            obj.result_E_dr(t_current)=  x(time * 7 + 1 + (IES_no - 1) * var );
            obj.result_H_dr(t_current) = x(time * 8 + 1 + (IES_no - 1) * var);
            %更新储能状态
            obj.ES_SOC(t_current+1) = obj.ES_selfd * obj.ES_SOC(t_current) - obj.result_ES_discharge(t_current) / obj.ES_eff / obj.ES_totalC + obj.result_ES_charge(t_current) * obj.ES_eff / obj.ES_totalC;
            obj.HS_SOC(t_current+1) = obj.HS_selfd * obj.HS_SOC(t_current) - obj.result_HS_discharge(t_current) / obj.HS_eff / obj.HS_totalC + obj.result_HS_charge(t_current) * obj.HS_eff / obj.HS_totalC;
        end
        
    end
    
end

