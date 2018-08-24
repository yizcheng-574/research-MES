classdef Grid_171118 < handle
    properties %可以设初始值
        %物理特性
        eleLimit_forward; % 大电网向微网群输电，正的
        eleLimit_reverse; % 回送，负的        
        %投标
        demand_curve;
        %最终优化结果
        result_demand;
    end
    
    methods
        function obj = Grid_171118(eleLimit_total) %初始化函数
            global priceNumbers period
            % 物理特性
            obj.eleLimit_forward = eleLimit_total(1); % 大电网向微网群输电
            obj.eleLimit_reverse = eleLimit_total(2); % 回送
            % 投标
            obj.demand_curve = zeros(priceNumbers, 1); %初始化投标曲线
            %最终优化结果
            obj.result_demand = zeros(24*period,1);
        end
        
            
        function gridDemand = zGenerate(obj, elePrice_current)
            global minMarketPrice priceNumbers step
            
            for i = 1 : priceNumbers
                pricePoint = minMarketPrice + (i-1) * step;
                
                if pricePoint <= elePrice_current
                    obj.demand_curve(i) = - obj.eleLimit_reverse; % 电价较低，电网购电，相当于本地负荷
                else
                    obj.demand_curve(i) = - obj.eleLimit_forward; % 电价较高，电网发电，相当于本地电源
                end
            end
            gridDemand = obj.demand_curve;
        end
        
        
        function clearDemand = getClearDemand(obj, clearPrice, t_current)
            global minMarketPrice step
            
            p = (clearPrice - minMarketPrice) / step + 1;
            p1 = floor((clearPrice - minMarketPrice) / step) + 1;
            p2 = ceil((clearPrice - minMarketPrice) / step) + 1;
            
            if p1 ~= p2
                slope = (obj.demand_curve(p2) - obj.demand_curve(p1)) / (p2 - p1);
                clearDemand = slope * (p - p1) + obj.demand_curve(p1);
            else
                clearDemand = obj.demand_curve(p1);
            end
            
            obj.result_demand(t_current) = clearDemand;
        end
        
        
        function gridDemand = handlePrice(obj, Eprice, t_current)
            global elePrice
            if Eprice <= elePrice(t_current)
                gridDemand = - obj.eleLimit_reverse; % 电价较低，电网购电，相当于本地负荷
            else
                gridDemand = - obj.eleLimit_forward; % 电价较高，电网发电，相当于本地电源
            end
        end
                
        
        
         % 输出优化结果
        function demandRecord = getResult(obj)
            demandRecord = obj.result_demand;
        end
        
        
    end
end