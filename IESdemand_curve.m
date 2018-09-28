function [demand,price] = IESdemand_curve(priceArray, pt)
    global EH1 EH2 EH3 IESNUMBER minMarketPrice gasPrice1
    price = minMarketPrice : 0.01 : 1;
    l = length(price);
    num = 1;
    demand = zeros( IESNUMBER, l);
    for p1 = minMarketPrice : 0.01 : 1
        priceArray(pt) = p1;
        for ies_no = 1: IESNUMBER
            eval(['[x,~,~,~,~] = EH',num2str(ies_no),'.handlePrice(priceArray, gasPrice1, pt);']);
            demand(ies_no,num) = x(1);
        end
        num = num + 1 ;
    end
%     figure
%     hold on
%     plot(price,demand1,'LineWidth',1.5);
%     plot(price,demand2,'LineWidth',1.5); 
%     plot(price,demand3,'LineWidth',1.5);
%     plot(price,demand1 + demand2 + demand3,'LineWidth',1.5);
%     xlabel('电价')
%     legend('IES1','IES2','IES3','总')
%     ylabel('需求')
end