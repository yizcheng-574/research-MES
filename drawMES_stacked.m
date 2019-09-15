function [c4_gridClearDemand] = drawMES_stacked(t, result_Ele,EH_res_total, ymin, ymax, limit)
    c4_gridClearDemand = sum(result_Ele , 2) - EH_res_total;   
    result_Ele_positive = result_Ele;
    result_Ele_positive(result_Ele_positive < 0) = 0; 
    result_Ele_negative = result_Ele;
    result_Ele_negative(result_Ele_negative > 0) = 0;
 
    figure;
    H1 = bar(t, result_Ele_positive / 1000, 'stacked'); hold on;
    H2 = bar(t, [result_Ele_negative, - EH_res_total]/1000, 'stacked');
    H1(1).EdgeColor = 'none';
    H1(2).EdgeColor = 'none';
    H1(3).EdgeColor = 'none';
    H2(1).EdgeColor = 'none';
    H2(2).EdgeColor = 'none';
    H2(3).EdgeColor = 'none';
    H2(4).EdgeColor = 'none';
    color_mes1 = ColorHex('4083ff') / 255;
    color_mes2 = ColorHex('005aff') / 255;
    color_mes3 = ColorHex('3200ff') / 255;
    gold = [1 0.843 0];

    H1(1).FaceColor = color_mes1;
    H1(2).FaceColor = color_mes2;
    H1(3).FaceColor = color_mes3;
    H2(1).FaceColor = color_mes1;
    H2(2).FaceColor = color_mes2;
    H2(3).FaceColor = color_mes3;
    H2(4).FaceColor = gold;
    P = plot(t, c4_gridClearDemand/1000);
    set(P, 'Color','Black', 'LineWidth', 1.5);
    if nargin > 5
       plot([0, t], ones(25, 1) * limit/1000, 'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',1);
    end
    ylabel('power(MW)');
    xlim([0, 25]);
    ylim([ymin, ymax]/1000)
    xticks(0: 6 : 24);
    xticklabels({ '0:00','6:00','12:00','18:00','24:00' });
    le = legend([H1(1), H1(2), H1(3), H2(4), P], ...
        'MES_1', 'MES_2', 'MES_3', 'RES integrated into the IMES', 'main transformer',...
        'Orientation','horizontal');
    set(le,'Box','off');
    set(gcf,'Position',[0 0 500 200]);

end

