function [] = drawMES(t, result_Ele, ymin, ymax)
  
    figure;
    
    H1 = bar(t, result_Ele / 1000); hold on;
    H1(1).EdgeColor = 'none';
    H1(2).EdgeColor = 'none';
    H1(3).EdgeColor = 'none';
   
    color_mes1 = ColorHex('4083ff') / 255;
    color_mes2 = ColorHex('005aff') / 255;
    color_mes3 = ColorHex('3200ff') / 255;

    H1(1).FaceColor = color_mes1;
    H1(2).FaceColor = color_mes2;
    H1(3).FaceColor = color_mes3;
   
    le = legend([H1(1), H1(2), H1(3)],'MES_1', 'MES_2', 'MES_3','Orientation','horizontal');
    ylabel('power(MW)');
    xlim([0, 25]);
    ylim([ymin, ymax]/1000);
    xticks(0: 6 : 24);
    xticklabels({ '0:00','6:00','12:00','18:00','24:00' });
    
    set(le,'Box','off');
    set(gcf,'Position',[0 0 500 200]);

end

