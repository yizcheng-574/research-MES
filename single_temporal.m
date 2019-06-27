%单时段的优化问题，与all_temporal比较
global minMarketPrice maxMarketPrice period IESNUMBER elePrice  minimumPower
ee = 1e-2; %0.0001 0.0003
iterativeStep = 1e-5; %0.00001 0.0001
iterationTimes = zeros(24*period, 2); %记录迭代次数
maxIteration = 3000; %最大迭代次数
gridClearDemand = zeros(24*period,1);
dflag = 0; % 是否绘图
for pt =  1 : 24 * period
     if isDA == 0
        EH1.predict(pt);
        EH2.predict(pt);
        EH3.predict(pt);
     end 
     if dflag == 1
        [demand,price] = IESdemand_curve(priceArray, pt);
        figure;
        hold on;
        for ies_no = 1 : IESNUMBER
            plot(price,demand(ies_no,:),'LineWidth',1.5);
        end
        plot(price,sum(demand),'LineWidth',1.5);
        xlabel('电价')
%             legend('IES1','IES2','IES3','总')
        title([pt ,'时刻投标曲线']);
        ylabel('需求')
     end

    % 最低价格，一般都是需求大于供给
    priceArray(pt) = minMarketPrice;
    [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, pt);
    clearDemand_minPrice_EH1 = x(1);
    [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, pt);
    clearDemand_minPrice_EH2 = x(1);
    [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice1, pt);
    clearDemand_minPrice_EH3 = x(1);
    if priceArray(pt) ==  elePrice(pt)
        clearDemand_minPrice_grid = clearDemand_minPrice_EH1 + clearDemand_minPrice_EH2 + clearDemand_minPrice_EH3;
        if clearDemand_minPrice_grid > eleLimit_total(1)
            clearDemand_minPrice_grid = eleLimit_total(1);
        end
        if clearDemand_minPrice_grid < minimumPower
            clearDemand_minPrice_grid = minimumPower;
        end
    elseif  priceArray(pt)>  elePrice(pt)
        clearDemand_minPrice_grid =eleLimit_total(1);
    else
        clearDemand_minPrice_grid =minimumPower;
    end
    clearDemand_minPrice = [-clearDemand_minPrice_grid; clearDemand_minPrice_EH1; clearDemand_minPrice_EH2; clearDemand_minPrice_EH3]; % 需求为正，供给为负

    % 最高价格，一般都是供给大于需求
    priceArray(pt) = maxMarketPrice;
    [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, pt);
    clearDemand_maxPrice_EH1 = x(1);
    [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, pt);
    clearDemand_maxPrice_EH2 = x(1);
    [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice1, pt);
    clearDemand_maxPrice_EH3 = x(1);
    if priceArray(pt) ==  elePrice(pt)
        clearDemand_maxPrice_grid = clearDemand_maxPrice_EH1 + clearDemand_maxPrice_EH2 + clearDemand_maxPrice_EH3;
        if clearDemand_maxPrice_grid > eleLimit_total(1)
            clearDemand_maxPrice_grid = eleLimit_total(1);
        end
        if clearDemand_maxPrice_grid < minimumPower
            clearDemand_maxPrice_grid = minimumPower;
        end
    elseif  priceArray(pt )>  elePrice(pt)
            clearDemand_maxPrice_grid =eleLimit_total(1);
    else
            clearDemand_maxPrice_grid =minimumPower;
    end
    clearDemand_maxPrice = [-clearDemand_maxPrice_grid; clearDemand_maxPrice_EH1; clearDemand_maxPrice_EH2; clearDemand_maxPrice_EH3]; % 需求为正，供给为负

    iterationNumber = 2;
    if sum(clearDemand_minPrice) * sum(clearDemand_maxPrice) <= 0 % 说明出清点在这个区间内，有两个问题，一是等于零是否直接结束，二是如果出清点不唯一怎么办
        % 市场出清得到出清价格，并更新预测电价序列
        [priceArray(pt), clearDemand] = iterativeClear(minMarketPrice, maxMarketPrice, clearDemand_minPrice, clearDemand_maxPrice, ee, priceArray, gasPrice1, pt);
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
[result_Ele(:,1), result_CHP_G(:,1), result_Boiler_G(:,1), result_ES_discharge(:,1), result_ES_charge(:,1), result_HS_discharge(:,1), result_HS_charge(:,1), result_ES_SOC(:,1), result_HS_SOC(:,1), result_EH_Le(:, 1), result_EH_Lh(:,1), result_EH_solarP(:,1), result_EH_windP(:,1), result_EH_Edr(:,1), result_EH_Hdr(:,1)] = EH1.getResult;
[result_Ele(:,2), result_CHP_G(:,2), result_Boiler_G(:,2), result_ES_discharge(:,2), result_ES_charge(:,2), result_HS_discharge(:,2), result_HS_charge(:,2), result_ES_SOC(:,2), result_HS_SOC(:,2), result_EH_Le(:, 2), result_EH_Lh(:,2), result_EH_solarP(:,2), result_EH_windP(:,2), result_EH_Edr(:,2), result_EH_Hdr(:,2)] = EH2.getResult;
[result_Ele(:,3), result_CHP_G(:,3), result_Boiler_G(:,3), result_ES_discharge(:,3), result_ES_charge(:,3), result_HS_discharge(:,3), result_HS_charge(:,3), result_ES_SOC(:,3), result_HS_SOC(:,3), result_EH_Le(:, 3), result_EH_Lh(:,3), result_EH_solarP(:,3), result_EH_windP(:,3), result_EH_Edr(:,3), result_EH_Hdr(:,3)] = EH3.getResult;

priceArray_record(:,3) = priceArray;