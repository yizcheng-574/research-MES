function [clearPrice, clearDemand] = iterativeClear(price_min, price_max, clearDemand_minPrice, clearDemand_maxPrice, ee, priceArray, gasPrice1, t_current)

    global Grid1 EH1 EH2 EH3 minMarketPrice maxMarketPrice iterationNumber elePrice eleLimit_total
    
    price_medium = (price_min + price_max)/2;
    
    priceArray(t_current) = price_medium;
    [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, t_current);
    clearDemand_medium_EH1 = x(1);
    [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, t_current);
    clearDemand_medium_EH2 = x(1);
    [x,~,~,~,~] = EH3.handlePrice(priceArray, gasPrice1, t_current);
    clearDemand_medium_EH3 = x(1);
    if priceArray(t_current) ==  elePrice(t_current)
        clearDemand_medium_grid = clearDemand_mediumPrice_EH1(i) + clearDemand_mediumPrice_EH2(i) + clearDemand_mediumPrice_EH3(i);
        if clearDemand_medium_grid > eleLimit_total(1)
            clearDemand_medium_grid = eleLimit_total(1);
        end
        if clearDemand_medium_grid < eleLimit_total(2)
            clearDemand_medium_grid = eleLimit_total(2);
        end
    elseif  priceArray(t_current)>  elePrice(t_current)
        clearDemand_medium_grid =eleLimit_total(1);
    else
        clearDemand_medium_grid =eleLimit_total(2);
    end
    clearDemand_mediumPrice = [-clearDemand_medium_grid; clearDemand_medium_EH1; clearDemand_medium_EH2; clearDemand_medium_EH3]; %需求为正，供给为负
    
    iterationNumber = iterationNumber + 1;
    
    if sum(clearDemand_mediumPrice) * sum(clearDemand_minPrice) <= 0 % 说明出清点在前半个区间
        distance = price_medium - price_min;
        if distance>ee
            % price_max = price_medium;
            [clearPrice, clearDemand] = iterativeClear(price_min, price_medium, clearDemand_minPrice, clearDemand_mediumPrice, ee, priceArray, gasPrice1,  t_current);
        else
            % 方法一，正常二分法的结果，但是不能保证供需平衡
            % clearPrice = (price_min + price_medium)/2; % 停止迭代
            
            % 方法二，假设线性，求出清价格和出清功率
            if sum(clearDemand_minPrice) == sum(clearDemand_mediumPrice)
                clearPrice = (price_min + price_medium)/2;
            else
                slope = (sum(clearDemand_minPrice) - sum(clearDemand_mediumPrice)) / (price_min - price_medium);
                clearPrice = price_min + (0 - sum(clearDemand_minPrice)) / slope;
            end
            
            if clearPrice < minMarketPrice || clearPrice > maxMarketPrice
                error('出清电价超出允许范围')
            end
            
            clearDemand = zeros(length(clearDemand_minPrice) ,1);
            for ii=1:length(clearDemand_minPrice)
                slope = (clearDemand_minPrice(ii) - clearDemand_mediumPrice(ii)) / (price_min - price_medium);
                clearDemand(ii) = clearDemand_minPrice(ii) + (clearPrice - price_min) * slope;
            end
            
        end
    else % 说明出清点在后半个区间
        distance = price_max - price_medium;
        if distance>ee
            % price_min = price_medium;
            [clearPrice, clearDemand] = iterativeClear(price_medium, price_max, clearDemand_mediumPrice, clearDemand_maxPrice, ee, priceArray, gasPrice1, t_current);
        else
            % 方法一，正常二分法的结果，但是不能保证供需平衡
            % clearPrice = (price_medium + price_max)/2; % 停止迭代            
            
            % 方法二，假设线性，求出清价格和出清功率
            if sum(clearDemand_mediumPrice) == sum(clearDemand_maxPrice)
                clearPrice = (price_medium + price_max)/2;
            else
                slope = (sum(clearDemand_mediumPrice) - sum(clearDemand_maxPrice)) / (price_medium - price_max);
                clearPrice = price_medium + (0 - sum(clearDemand_mediumPrice)) / slope;
            end
            
            if clearPrice < minMarketPrice || clearPrice > maxMarketPrice
                error('出清电价超出允许范围')
            end
            
            clearDemand = zeros(length(clearDemand_mediumPrice) ,1);
            for ii=1:length(clearDemand_mediumPrice)
                slope = (clearDemand_mediumPrice(ii) - clearDemand_maxPrice(ii)) / (price_medium - price_max);
                clearDemand(ii) = clearDemand_mediumPrice(ii) + (clearPrice - price_medium) * slope;
            end
            
        end
    end
   
end