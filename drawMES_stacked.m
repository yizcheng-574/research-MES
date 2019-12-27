function [c4_gridClearDemand] = drawMES_stacked(t, result_Ele, EH_res_total, ymin, ymax, datafile, elePrice)
    global isEn
    load(datafile, 'eleLimit_total', 'priceArray');
    plotAux;    
    c4_gridClearDemand = sum(result_Ele , 2) - EH_res_total;   
    result_Ele_positive = result_Ele;
    result_Ele_positive(result_Ele_positive < 0) = 0; 
    result_Ele_negative = result_Ele;
    result_Ele_negative(result_Ele_negative > 0) = 0;
 
    figure;
    yyaxis left;
    H1 = bar(t, result_Ele_positive / 1000, 'stacked'); hold on;
    H2 = bar(t, [result_Ele_negative, - EH_res_total]/1000, 'stacked');
    H1(1).EdgeColor = 'none';
    H1(2).EdgeColor = 'none';
    H1(3).EdgeColor = 'none';
    H2(1).EdgeColor = 'none';
    H2(2).EdgeColor = 'none';
    H2(3).EdgeColor = 'none';
    H2(4).EdgeColor = 'none';
    
    color_mes1 = [0.5, 0.5, 0.5];
    color_mes2 = ColorHex('4083ff') / 255;
%     color_mes2 = ColorHex('005aff') / 255;
    color_mes3 = ColorHex('3200ff') / 255;

    H1(1).FaceColor = color_mes1;
    H1(2).FaceColor = color_mes2;
    H1(3).FaceColor = color_mes3;
    H2(1).FaceColor = color_mes1;
    H2(2).FaceColor = color_mes2;
    H2(3).FaceColor = color_mes3;
    H2(4).FaceColor = gold;
    
    P = plot(t, c4_gridClearDemand/1000);
    set(P, 'Color', 'Black', 'LineWidth', 1.5);
    
    plot([0, t], ones(25, 1) * eleLimit_total(1)/1000, 'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',1);
    plot([0, t], ones(25, 1) * eleLimit_total(2)/1000, 'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',1);

    if isEn == 1
        ylabel('power(MW)');
    else
        ylabel('功率(MW)');
    end
    ylim([ymin, ymax]/1000);
    
    yyaxis right;
    
    H3 = plot(1:24, [priceArray, elePrice]);
    set(H3(1),'Color',firebrick, 'LineStyle','--','LineWidth',1.5);
    set(H3(2),'Color',firebrick, 'LineStyle','-','LineWidth',1.5);
    if isEn == 1
        ylabel('electricity price(yuan/kWh)');
    else
        ylabel('电价(元/kWh)');
    end
    ylim([0.1, 1.2]);yticks(0.1 : 0.3 : 1.2);
    
    xlim([0, 25]); xticks(0: 6 : 24); xticklabels({ '0:00','6:00','12:00','18:00','24:00' });
    if isEn == 1
        le = legend([H1(1), H1(2), H1(3), H2(4), P, H3(1), H3(2)], ...
            'MES_1', 'MES_2', 'MES_3', 'shared RES', 'main transformer',...
            'clearing price','RTP price',...
            'Orientation','horizontal');
    else
        le = legend([H1(1), H1(2), H1(3), H2(4), P, H3(1), H3(2)], ...
            'MES_1', 'MES_2', 'MES_3', 'IMES共享可再生能源', '主变压器',...
            '出清电价','主网实时电价',...
            'Orientation','horizontal');
         xlabel('时间')
    end
    set(le,'Box','off');
%     set(le, 'NumColumns', 4);
    set(gcf,'Position',[0 0 500 250]);

end