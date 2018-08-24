% clc
% clear

%多能协调优化
%能源枢纽的本地集中式优化策略
%2017.8.22 ver1.0
%2017.8.28 ver1.1 电、热功率平衡改为大于约束，增加购气量和的约束，BUG为解决：热有盈余的话储热会又充又放
%2017.8.28 ver2.0 改为class
%2018.1.10 ver3.0 market信息改用全局变量，对储能的又充又放采用等效的线性优化方法，时间间隔改为15min

% f2=[-13;-23];
% A2=[5,15; 4,4; 35,20;];
% b2=[480,160,1190];
% lb=[0;0];
% [x,fval,exitflag,output,lambda] = linprog(f2,A2,b2,[],[],lb,[])



classdef EH_local_170828_v3 < handle
    properties %可以设初始值
        %负荷与可再生能源
        Le;
        Lh;
        solarP;
        solarP_rate;
        windP;
        windP_rate;
        dev_L; %百分数
        dev_PV;
        dev_WT;
        
        %CHP的参数
        CHP_GE_eff;
        CHP_GH_eff;
        CHP_G_max;
        CHP_G_min;
        
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
        
        %热储能
        HS_totalC;
        HS_maxSOC;
        HS_minSOC;
        % HS_currentC; %随时间变化的
        HS_targetSOC;
        HS_Hmax;
        HS_eff;
        
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
    end
    
    methods
        
        function obj = EH_local_170828_v3(eleLimit, gasLimit, Le_base, Lh_base, solar_base, wind_base, CHP_para, Boiler_para, ES_para, HS_para, dev_load, dev_solar, dev_wind, solar_rate, wind_rate)
            
            global priceNumbers period
            
            %联络线电交换和气交换约束
            obj.Ele_max = eleLimit(1);
            obj.Ele_min = eleLimit(2);
            obj.Ele_eff = eleLimit(3); %购电效率，即线路传输的线损率
            obj.Gas_max = gasLimit; %两个购气量的和
            % Gas_min = 0;
            
            %负荷与可再生能源
            obj.Le = Le_base; %base是日前初步的预测曲线
            obj.Lh = Lh_base;
            obj.solarP = solar_base;
            obj.solarP_rate = solar_rate;
            obj.windP = wind_base;
            obj.windP_rate = wind_rate;
            obj.dev_L = dev_load; %百分数
            obj.dev_PV = dev_solar;
            obj.dev_WT = dev_wind;
            
            %CHP的参数
            obj.CHP_GE_eff = CHP_para(1);
            obj.CHP_GH_eff = CHP_para(2);
            obj.CHP_G_max = CHP_para(3) / obj.CHP_GE_eff; %由额定电功率得到最大购气量
            obj.CHP_G_min = CHP_para(4) * obj.CHP_G_max;
            
            %锅炉
            obj.Boiler_eff = Boiler_para(1);
            obj.Boiler_G_max = Boiler_para(2) / obj.Boiler_eff;  %10
            % Boiler_G_min = 0;
            
            %电储能和热储能
            obj.ES_totalC = ES_para(1);
            obj.ES_maxSOC = ES_para(2);
            obj.ES_minSOC = ES_para(3);
            obj.ES_targetSOC = ES_para(5);
            obj.ES_Pmax = obj.ES_totalC / ES_para(6);
            obj.ES_eff = ES_para(7);
            
            obj.HS_totalC = HS_para(1);
            obj.HS_maxSOC = HS_para(2);
            obj.HS_minSOC = HS_para(3);
            obj.HS_targetSOC = HS_para(5);
            obj.HS_Hmax = obj.HS_totalC / HS_para(6);
            obj.HS_eff = HS_para(7);
            
            % 投标
            obj.demand_curve = zeros(priceNumbers, 1); %初始化投标曲线
            
            %最终优化结果
            obj.ES_SOC = zeros(24*period+1,1);
            obj.ES_SOC(1) = ES_para(4);
            obj.HS_SOC = zeros(24*period+1,1);
            obj.HS_SOC(1) = HS_para(4);
            obj.result_Ele = zeros(24*period,1);
            obj.result_CHP_G = zeros(24*period,1);
            obj.result_Boiler_G = zeros(24*period,1);
            obj.result_ES_discharge = zeros(24*period,1);
            obj.result_ES_charge = zeros(24*period,1);
            obj.result_HS_discharge = zeros(24*period,1);
            obj.result_HS_charge = zeros(24*period,1);
        end
        
        % 可再生能源与负荷的预测业务，一天预测25次，日前一次，日内每小时一次
        function predict(obj, t_current) % 第几小时的预测，time=0是日前，1-24是日内
            %             global period
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
                Le_error = randn() * obj.Le(t_current) * obj.dev_L; %预测误差
                Lh_error = randn() * obj.Lh(t_current) * obj.dev_L;
                obj.Le(t_current) = obj.Le(t_current) + Le_error; %更新
                obj.Lh(t_current) = obj.Lh(t_current) + Lh_error;
                
                %风、光
                %有几点和负荷不一样
                %一是：如果光是零，那么一定是零，没有预测误差
                %二是：风、光接近零的时候，不要随机变为负数
                %                 solarP_error = zeros(24*period,1);
                %                 windP_error = zeros(24*period,1);
                %                 solarP_error(t_current:24*period) = randn([(24*period-t_current+1),1]) .* obj.solarP(t_current:24*period) * obj.dev_RES; %预测误差
                %                 windP_error(t_current:24*period) = randn([(24*period-t_current+1),1]) .* obj.windP(t_current:24*period) * obj.dev_RES;
                %                 %修正
                %                 for i = t_current : 24*period
                %                     if obj.solarP(i) == 0
                %                         solarP_error(i) = 0;
                %                     end
                %                     if obj.solarP(i) + solarP_error(i) < 0
                %                         solarP_error(i) = - obj.solarP(i);
                %                     end
                %                     if obj.windP(i) == 0
                %                         windP_error(i) = 0;
                %                     end
                %                     if obj.windP(i) + windP_error(i) < 0
                %                         windP_error(i) = - obj.windP(i);
                %                     end
                %                 end
                %                 obj.solarP = obj.solarP + solarP_error; %更新
                %                 obj.windP = obj.windP + windP_error;
                %改为只计当前时刻的预测误差，否则累计的误差过大；频繁重复预测也不科学
                solarP_error = randn() * obj.solarP(t_current) * obj.dev_PV; %预测误差，solarP(t_current)=0的时候自然为零
                windP_error = randn() * obj.windP(t_current) * obj.dev_WT; %预测误差
                
                obj.solarP(t_current) = obj.solarP(t_current) + solarP_error; %更新
                if obj.solarP(t_current) < 0
                    obj.solarP(t_current) = 0;
                elseif obj.solarP(t_current) > obj.solarP_rate
                    obj.solarP(t_current) = obj.solarP_rate;
                end
                
                obj.windP(t_current) = obj.windP(t_current) + windP_error;
                if obj.windP(t_current) < 0
                    obj.windP(t_current) = 0;
                elseif obj.windP(t_current) > obj.windP_rate
                    obj.windP(t_current) = obj.windP_rate;
                end
            else
                %日前预测值就是base值
            end
        end
        
        
        
        %通过多次解最优化，构造投标函数
        function demand_curve_result = curveGenerate(obj, Eprice, Gprice, t_current) %接收电价、气价、当前时间
            % predict(obj, t_current); %更新本地的负荷预测
            global minMarketPrice step priceNumbers
            
            for i = 1: 1 : priceNumbers
                Eprice(t_current) = minMarketPrice + (i-1) * step;
                [x,~,~,~,~] = localOptimal(obj, Eprice, Gprice, t_current, 9e9); % conditionEle = 9e9
                obj.demand_curve(i) = x(1);
            end
            
            demand_curve_result = obj.demand_curve;
        end
        
        
        %接收市场出清价格，做自治优化，更新自身状态
        function [x,fval,exitflag,output,lambda] = handlePrice(obj, Eprice, Gprice, t_current) %这里的x随着t_current会越来越少
            global period
            
            [x,fval,exitflag,output,lambda] = localOptimal(obj, Eprice, Gprice, t_current, 9e9); % conditionEle = 9e9
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
        
        
        % 用于TC单次出清
        %接收市场出清价格，按投标曲线得到出清功率，在保持出清功率不变的情况下，求解最优化
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
        
        
        % 用于TC迭代出清 二分法
        % 接收市场出清价格、出清功率，在保持出清功率不变的情况下，求解最优化
        function [x,fval,exitflag,output,lambda] = conditionHandlePrice_2(obj, Eprice, Gprice, t_current, clearDemand) %这里的x随着t_current会越来越少
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
            
            %更新储能状态
            obj.ES_SOC(t_current+1) = obj.ES_SOC(t_current) - obj.result_ES_discharge(t_current) / obj.ES_eff / obj.ES_totalC + obj.result_ES_charge(t_current) * obj.ES_eff / obj.ES_totalC;
            obj.HS_SOC(t_current+1) = obj.HS_SOC(t_current) - obj.result_HS_discharge(t_current) / obj.HS_eff / obj.HS_totalC + obj.result_HS_charge(t_current) * obj.HS_eff / obj.HS_totalC;
        end
        
        
        
        function [x,fval,exitflag,output,lambda] = conditionHandlePrice_DA(obj, Eprice, Gprice, t_current, clearDemand) %这里的x随着t_current会越来越少
            global period
            conditionEle = clearDemand;
            
            [x,fval,exitflag,output,lambda] = localOptimal(obj, Eprice, Gprice, t_current, conditionEle);
            time = 24*period - t_current + 1; %总时间段
            for pt = t_current: 24 * period
            % 只执行当前周期的结果
            obj.result_Ele(pt) = x(1 + pt - t_current);
            obj.result_CHP_G(pt) = x(time + 1 + pt - t_current);
            obj.result_Boiler_G(pt) = x(time*2 + 1 + pt - t_current);
            obj.result_ES_discharge(pt) = x(time*3 + 1 + pt - t_current);
            obj.result_ES_charge(pt) = x(time*4 + 1 + pt - t_current);
            obj.result_HS_discharge(pt) = x(time*5 + 1 + pt - t_current);
            obj.result_HS_charge(pt) = x(time*6 + 1 + pt - t_current);
            
            %更新储能状态
            obj.ES_SOC(pt+1) = obj.ES_SOC(pt) - obj.result_ES_discharge(pt) / obj.ES_eff / obj.ES_totalC + obj.result_ES_charge(pt) * obj.ES_eff / obj.ES_totalC;
            obj.HS_SOC(pt+1) = obj.HS_SOC(pt) - obj.result_HS_discharge(pt) / obj.HS_eff / obj.HS_totalC + obj.result_HS_charge(pt) * obj.HS_eff / obj.HS_totalC;
        
            end
        end
        % 输出优化结果
        function [Ele, G_CHP, G_Boiler, ES_discharge, ES_charge, HS_discharge, HS_charge, ES_SOC, HS_SOC, Le, Lh, solarP, windP] = getResult(obj)
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
            Le = obj.Le;
            Lh = obj.Lh;
            solarP = obj.solarP;
            windP = obj.windP;
        end
        
        
        
        % 测试优化结果
        function [result_balance_P, result_balance_H, result_check_ES, result_check_HS] = testResult(obj)
            %电、热功率平衡性测试，降低要求至大于零
            result_balance_P = obj.result_Ele + obj.CHP_GE_eff.*obj.result_CHP_G + obj.result_ES_discharge - obj.result_ES_charge - obj.Le + obj.windP + obj.solarP;
            result_balance_H = obj.CHP_GH_eff.*obj.result_CHP_G + obj.Boiler_eff.*obj.result_Boiler_G + obj.result_HS_discharge - obj.result_HS_charge - obj.Lh;
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
            time = 24*period - t_current + 1; %总时间段
            var = time * 7; %总变量数
            %第1,2,3组time是购电量、CHP购气量、锅炉购气量，第4-7组time是储电、储热的放、充功率
            
            
            
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
                ub(time+i, 1) = obj.CHP_G_max;
                ub(time*2+i, 1) = obj.Boiler_G_max;
                ub(time*3+i, 1) = obj.ES_Pmax;
                ub(time*4+i, 1) = obj.ES_Pmax;
                ub(time*5+i, 1) = obj.HS_Hmax;
                ub(time*6+i, 1) = obj.HS_Hmax;
            end
            for i = 1 : time
                lb(i, 1) = obj.Ele_min;
                lb(time+i, 1) = obj.CHP_G_min;
                %                     lb(time*2+i, 1) = 0;
                %                     lb(time*3+i, 1) = 0;
                %                     lb(time*4+i, 1) = 0;
                %                     lb(time*5+i, 1) = 0;
                %                     lb(time*6+i, 1) = 0;
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
            %beq_Ebus = - obj.Le;
            %beq_Hbus = - obj.Lh;
            beq_Ebus = - obj.Le(t_current : 24*period) + obj.windP(t_current : 24*period) + obj.solarP(t_current : 24*period); %Le的size不变，但是取部分值
            beq_Hbus = - obj.Lh(t_current : 24*period); %Lh的size不变，但是取部分值
            for i=1:time
                Aeq_Ebus(i,i) = - obj.Ele_eff; %线损率
                Aeq_Ebus(i,time+i) = - obj.CHP_GE_eff;
                Aeq_Ebus(i,time*3+i) = - 1; %放电
                Aeq_Ebus(i,time*4+i) = 1; %充电
            end
            for i=1:time
                Aeq_Hbus(i,time+i) = - obj.CHP_GH_eff;
                Aeq_Hbus(i,time*2+i) = - obj.Boiler_eff;
                Aeq_Hbus(i,time*5+i) = - 1; %放热
                Aeq_Hbus(i,time*6+i) = 1; %充热
            end
            
            
            %电、热储能平衡性约束
            Aeq_ES = zeros(1, var);
            beq_ES = - (obj.ES_targetSOC - obj.ES_SOC(t_current)) * obj.ES_totalC;
            for i=1:time
                Aeq_ES(1, time*3+i) = 1/obj.ES_eff; %放电
                Aeq_ES(1, time*4+i) = - 1*obj.ES_eff; %充电
            end
            Aeq_HS = zeros(1, var);
            beq_HS = - (obj.HS_targetSOC - obj.HS_SOC(t_current)) * obj.HS_totalC;
            for i=1:time
                Aeq_HS(1, time*5+i) = 1/obj.HS_eff; %放热
                Aeq_HS(1, time*6+i) = - 1*obj.HS_eff; %充热
            end
            
            
            
            %不等式约束包括：SOC约束；购气量和的约束；（爬坡率约束）
            %SOC约束 A1是上限，A2是下限
            A1_Esoc = zeros(time, var);
            A2_Esoc = zeros(time, var);
            b1_Esoc = ones(time,1) * (obj.ES_maxSOC - obj.ES_SOC(t_current)) * obj.ES_totalC;
            b2_Esoc = ones(time,1) * (obj.ES_SOC(t_current) - obj.ES_minSOC) * obj.ES_totalC;
            for i=1:time
                for j=1 : i
                    A1_Esoc(i, time*3+j) = -1/obj.ES_eff; %放电
                    A1_Esoc(i, time*4+j) = 1*obj.ES_eff; %充电
                end
            end
            for i=1:time
                for j=1 : i
                    A2_Esoc(i, time*3+j) = 1/obj.ES_eff; %放电
                    A2_Esoc(i, time*4+j) = -1*obj.ES_eff; %充电
                end
            end
            
            A1_Hsoc = zeros(time, var);
            A2_Hsoc = zeros(time, var);
            b1_Hsoc = ones(time,1) * (obj.HS_maxSOC - obj.HS_SOC(t_current)) * obj.HS_totalC;
            b2_Hsoc = ones(time,1) * (obj.HS_SOC(t_current) - obj.HS_minSOC) * obj.HS_totalC;
            for i=1:time
                for j=1 : i
                    A1_Hsoc(i, time*5+j) = -1/obj.HS_eff; %放热
                    A1_Hsoc(i, time*6+j) = 1*obj.HS_eff; %充热
                end
            end
            for i=1:time
                for j=1 : i
                    A2_Hsoc(i, time*5+j) = 1/obj.HS_eff; %放热
                    A2_Hsoc(i, time*6+j) = -1*obj.HS_eff; %充热
                end
            end
            
            %购气量和的约束
            A_Gmax = zeros(time, var);
            b_Gmax = ones(time,1) .* obj.Gas_max;
            for i=1:time
                A_Gmax(i, time+i) = 1;
                A_Gmax(i, time*2+i) = 1;
            end
            
            
            
            %归纳所有线性约束
            %等式约束包括：电、热平衡约束（改为不等式），电、热储能平衡性约束（改为不等式）
            %不等式约束包括：SOC约束，购气量和的约束，（爬坡率约束）
            Aeq=[];
            beq=[];
            A=[Aeq_Ebus; Aeq_Hbus;   Aeq_ES; Aeq_HS;    A1_Esoc; A2_Esoc; A1_Hsoc; A2_Hsoc;  A_Gmax];
            b=[beq_Ebus; beq_Hbus;   beq_ES; beq_HS;    b1_Esoc; b2_Esoc; b1_Hsoc; b2_Hsoc;  b_Gmax];
            
            %             %fmincon需要列出一个初始可行解
            %             x0 = zeros(var,1);
            %             for i = 1 : time
            %                 x0(i) = min(-beq_Ebus(i), obj.Ele_max); % 购电量，优先线路购电
            %                 x0(time+i) = (-beq_Ebus(i) - x0(i)) / obj.CHP_GE_eff; % CHP购气量，不够的电由CHP发，顺便发热
            %                 x0(time*2+i) = max(-beq_Hbus(i) - x0(time+i)*obj.CHP_GH_eff , 0) / obj.Boiler_eff; % boiler购气量，不够的热由锅炉发
            %             %     x0(time*3+i, 1) = 0;
            %             %     x0(time*4+i, 1) = 0;
            %             %     x0(time*5+i, 1) = 0;
            %             %     x0(time*6+i, 1) = 0;
            %             end
            
            
            [x,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub);
            
            %             options = optimoptions('fmincon','MaxFunEvals',1000000);
            %             [x,fval,exitflag,output,lambda] = fmincon('myfun_1', x0, A, b, Aeq, beq, lb, ub, 'mycon', options);
            
            if exitflag ~= 1
                %                 error('没有可行解')
            end
            
        end
        
    end
    
end

