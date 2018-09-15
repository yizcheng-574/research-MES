%单时段的优化问题，与all_temporal比较
global minMarketPrice maxMarketPrice period
ee = 1e-2; %0.0001 0.0003
iterativeStep = 5e-6; %0.00001 0.0001
iterationTimes = zeros(24*period, 2); %记录迭代次数
maxIteration = 3000; %最大迭代次数
gridClearDemand = zeros(24*period,1);
if isDA 
   EH1.predict(0);
   EH2.predict(0);
   EH3.predict(0);
   priceArray_record(:,1) = elePrice;
end
if isGrad == 1%次梯度法求解
    for pt =  1 : 24 * period
        if isDA ~= 0
            EH1.predict(pt);
            EH2.predict(pt);
            EH3.predict(pt);
        end
         % 最低价格，一般都是需求大于供给
        priceArray(pt) = minMarketPrice;
        [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, pt);
        clearDemand_minPrice_EH1 = x(1);
        [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, pt);
        clearDemand_minPrice_EH2 = x(1);
        [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice1, pt);
        clearDemand_minPrice_EH3 = x(1);
        clearDemand_minPrice_grid = Grid1.handlePrice(priceArray(pt), pt);
        clearDemand_minPrice = [clearDemand_minPrice_grid; clearDemand_minPrice_EH1; clearDemand_minPrice_EH2; clearDemand_minPrice_EH3]; % 需求为正，供给为负
        
        % 最高价格，一般都是供给大于需求
        priceArray(pt) = maxMarketPrice;
        [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, pt);
        clearDemand_maxPrice_EH1 = x(1);
        [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, pt);
        clearDemand_maxPrice_EH2 = x(1);
        [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice1, pt);
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
                if isDA
                    break;
                else
                    error('超出最大迭代次数');
                end
            end
            
            %当前价格下的出力
            priceArray(pt) = elePrice(pt) + lamda_new;
            [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, pt);
            clearDemand_EH1_new = x(1);
            [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, pt);
            clearDemand_EH2_new = x(1);
            [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice1, pt);
            clearDemand_EH3_new = x(1);
            
            % clearDemand_grid_new = Grid1.handlePrice(priceArray(pt), pt);
            
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
            % lamda_new = max(0, lamda_old + sum(clearDemand) * iterativeStep);
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
            
            % lamdaArrayAvg_new = sum(lamdaArray) / length(lamdaArray); %方法一
            % lamda_avg_new = 1 / length(lamda_record) * lamda_record(number) + (length(lamda_record)-1) / length(lamda_record) * lamda_avg_old; %方法二
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
        EH3.conditionHandlePrice_2(priceArray, gasPrice1, pt, clearDemand(4));
        
    end
    
else    %二分法求解
    for pt =  1 : 24*period
        % 最低价格，一般都是需求大于供给
        priceArray(pt) = minMarketPrice;
        [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, pt);
        clearDemand_minPrice_EH1 = x(1);
        [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, pt);
        clearDemand_minPrice_EH2 = x(1);
        [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice1, pt);
        clearDemand_minPrice_EH3 = x(1);
        clearDemand_minPrice_grid = Grid1.handlePrice(priceArray(pt), pt);
        clearDemand_minPrice = [clearDemand_minPrice_grid; clearDemand_minPrice_EH1; clearDemand_minPrice_EH2; clearDemand_minPrice_EH3]; % 需求为正，供给为负
        
        % 最高价格，一般都是供给大于需求
        priceArray(pt) = maxMarketPrice;
        [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, pt);
        clearDemand_maxPrice_EH1 = x(1);
        [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, pt);
        clearDemand_maxPrice_EH2 = x(1);
        [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice1, pt);
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
        EH3.conditionHandlePrice_2(priceArray, gasPrice1, pt, clearDemand(4));
        
    end
end
[result_Ele(:,1), result_CHP_G(:,1), result_Boiler_G(:,1), result_ES_discharge(:,1), result_ES_charge(:,1), result_HS_discharge(:,1), result_HS_charge(:,1), result_ES_SOC(:,1), result_HS_SOC(:,1), EH1_Le, EH1_Lh, EH1_solarP, EH1_windP, EH1_Edr, EH1_Hdr] = EH1.getResult;
[result_Ele(:,2), result_CHP_G(:,2), result_Boiler_G(:,2), result_ES_discharge(:,2), result_ES_charge(:,2), result_HS_discharge(:,2), result_HS_charge(:,2), result_ES_SOC(:,2), result_HS_SOC(:,2), EH2_Le, EH2_Lh, EH2_solarP, EH2_windP, EH2_Edr, EH2_Hdr] = EH2.getResult;
[result_Ele(:,3), result_CHP_G(:,3), result_Boiler_G(:,3), result_ES_discharge(:,3), result_ES_charge(:,3), result_HS_discharge(:,3), result_HS_charge(:,3), result_ES_SOC(:,3), result_HS_SOC(:,3), EH3_Le, EH3_Lh, EH3_solarP, EH3_windP, EH3_Edr, EH3_Hdr] = EH3.getResult;
if isDA 
    priceArray_record(:,2) = priceArray;
else
    priceArray_record(:,3) = priceArray;
end