close all
% 数据处理
global period caseType
%--------------------颜色尝试---------------------
    orange = [1 0.65 0];
    gold = [1 0.843 0];
    gray = [0.5 0.5 0.5];
    
    olivedrab = [0.41961 0.55686 0.13725];
    yellowgreen = [0.60392 0.80392 0.19608];
    
    firebrick = [0.69804 0.13333 0.13333];
    tomato = [1 0.38824 0.27843];
    brown = [0.80392 0.2 0.2];
    maroon = [0.6902 0.18824 0.37647];
    
    royalblue = [0.2549 0.41176 0.88235];
    royalblue_dark = [0.15294 0.25098 0.5451];
    darkblue =[0 0 0.5451];
    dodgerblue = [0.11765 0.56471 1];
    
    indianred = [1 0.41 0.42];
    chocolate3 = [0.804 0.4 0.113];
    tan2 = [0.93  0.60 0.286];
    
    c1 = ColorHex('0D56A6');
    c2 = ColorHex('41DB00');
    c3 = ColorHex('A63C00');
%--------------------------------------负荷和可再生能源的曲线--------------------------------------
if caseType ~=32
    figure(1)
    optNumber=24;
    t=1:1:24*period;
    w=1.2;
    
    for IES_no = 1 : 3
        t_1 = 0 : 24;
        EH_Le_base = result_EH_Le(:, IES_no);

        EH_Lh_base = result_EH_Lh(:, IES_no);
        EH_solarP = result_EH_solarP(:, IES_no);
        EH_windP = result_EH_windP(:, IES_no);
        EH_Le = EH_Le_base + result_EH_Edr(:, IES_no);
        EH_Lh = EH_Lh_base + result_EH_Hdr(:, IES_no);
        subplot(3 , 2 , (IES_no - 1) * 2 + 1 )
        hold on;
        stairs(t_1, appendStairArray(EH_Le_base) / 1000, 'Color', 'b', 'LineStyle', '-.', 'LineWidth', w)
        stairs(t_1, appendStairArray(EH_Le) / 1000, 'Color', 'b', 'LineStyle', '-', 'LineWidth', w) 
        stairs(t_1, appendStairArray(EH_Lh_base) / 1000, 'Color', 'r', 'LineStyle', '-.', 'LineWidth',w)
        stairs(t_1, appendStairArray(EH_Lh) / 1000, 'Color', 'r', 'LineStyle', '-', 'LineWidth',w)
        xlim([0, 24 * period])
        ylim([0, max(max(EH_Le), max(EH_Lh)) / 1000]);
        xticks(0 : (24 * period / 4) : 24 * period);
        xticklabels({ '0:00', '6:00', '12:00', '18:00', '24:00' });
        xlabel(sprintf('MES%d', IES_no))
        ylabel('Load(MW)')
        if IES_no == 1
            H1 = legend('base electric load','total electric load','base thermal load','total thermal load',...
                'Location','northoutside','Orientation','horizontal');
            set(H1,'Box','off');
        end
        % ylabel('load / kW')
        % xlabel('time / h')
        % legend('Le','Lh','Location','northoutside','Orientation','horizontal')
        % xlabel('时间(h)')
        subplot(3,2,(IES_no - 1) * 2 + 2)
        hold on
        stairs(t_1, appendStairArray(EH_solarP) / 1000, 'Color', 'k', 'LineStyle', '-', 'LineWidth', w);
        stairs(t_1, appendStairArray(EH_windP) / 1000, 'Color', 'k', 'LineStyle', '--', 'LineWidth', w);
        xlim([0, 24*period]);
        xticks(0 : (24 * period /4) : 24 * period);
        xticklabels({'0:00','6:00','12:00','18:00','24:00'});
        xlabel(sprintf('MES%d', IES_no))
        % ylabel('RES power / kW')
        % xlabel('time / h')
        % legend('PV','WT','Location','northoutside','Orientation','horizontal')
        ylabel('power(MW)')
        if IES_no == 1
            le = legend('PV','wind','Location','northoutside','Orientation','horizontal');
            set(le,'Box','off');
        end
    end
    set(gcf,'Position',[0 0 800 500]);
    
    
    %--------------------------------------数据处理--------------------------------------
    result_Gas = result_CHP_G + result_Boiler_G;
    for IES_no = 1 : 3
        eval(['result_Ele_loss(:,IES_no) = result_Ele(:,IES_no) .* EH',num2str(IES_no),'.Ele_eff;']); % eleLimit(3)是线损率
        eval(['result_CHP_power(:,IES_no) = result_CHP_G(:,IES_no) .* EH',num2str(IES_no),'.CHP_GE_eff; ']);
        eval(['result_CHP_heat(:,IES_no) = result_CHP_G(:,IES_no) .* EH',num2str(IES_no),'.CHP_GH_eff; ']);
        eval(['result_Boiler_heat(:,IES_no) = result_Boiler_G(:,IES_no) .* EH',num2str(IES_no),'.Boiler_eff;']);
        eval(['result_eBoiler_H(:, IES_no) = result_eBoiler_E(:, IES_no) .* EH',num2str(IES_no),'.eBoiler_eff;']);
    end
    %--------------------------------------测试优化结果--------------------------------------
    ee = 1e-3;
    
    % 当储能的最大充、放电功率很大时，1000 * 1e-3 也会越线，因此应该提前将1e-3置为0
    result_ES_discharge(result_ES_discharge < ee) = 0;
    result_ES_charge(result_ES_charge < ee) = 0;
    result_HS_discharge(result_HS_discharge < ee) = 0;
    result_HS_charge(result_HS_charge < ee) = 0;
    
    % 计算
    %计算总成本 按网价计算
    if exist('priceArray', 'var')
        cost_clear =  (result_Ele' * priceArray + sum(result_Gas)' * gasPrice1) / period;
    end
    cost_utility = (result_Ele' *  elePrice + sum(result_Gas)' * gasPrice1) / period;
    
    totalCost = sum(cost_utility)
    % --------------------------------------绘图--------------------------------------
    t1 = 1 : 1 : 24 * period;
    t2 = 0 : 1 : 24 * period;
    optNumber = 24;
    w=1.5;
    
    %--------------------阻塞管理---------------------
    if isCentral == 0
        c4_clearingPrice = priceArray;
        c4_gridClearDemand = - sum(result_Ele , 2);
    else
        c4_clearingPrice = elePrice;
        c4_gridClearDemand = - sum(result_Ele , 2);
    end
    figure;
    hold on;
    yyaxis left;
    H1 = bar(t1, -c4_gridClearDemand/1000);
    H1(1).FaceColor = dodgerblue;
    H1(1).EdgeColor = 'none';
    H3 = stairs(t2, ones(24*period+1, 1) .* eleLimit_total(1)/1000, 'Color',gray,'LineStyle','--','LineWidth',1);
    H4 = stairs(t2, ones(24*period+1, 1) * eleLimit_total(2)/1000, 'Color',gray,'LineStyle','--','LineWidth',1);
    ylabel('transformer power(MW)');
    uplimit= ceil(eleLimit_total(1) / 1000 * 1.1);
    lowerlimit=ceil(-eleLimit_total(2) / 1000 * 1.1);
    ylim([-lowerlimit, uplimit]);
    yticks(-lowerlimit : 1 : uplimit);
    
    yyaxis right;
    H2 = plot(t1, [c4_clearingPrice, elePrice]);
    le = legend([H1,H2(1),H2(2),H3(1)],...
        'transformer power','hourly clearing price','utility price','power limits');...
    set(le,'Box','off');
    set(H2(1),'Color',firebrick, 'LineStyle','-','LineWidth',1.5, 'Marker', '.', 'MarkerSize', 13)
    set(H2(2),'Color',darkblue, 'LineStyle','-','LineWidth',1.5, 'Marker', '.', 'MarkerSize', 13)
    ylabel('electricity price(yuan/kWh)');
    ylim([0.1, 1.2]);
    yticks(0.1 : 0.3 : 1.2);
    %xlabel('时间(h)')
    xlim([0, 24 * period + 1]);
    xticks(0 : (24*period/4):24*period);
    xticklabels({'0:00','6:00','12:00','18:00','24:00'});
    set(gcf,'Position',[0 0 400 200]);
    
    %--------------------优化结果2---------------------
    result_Ele_loss_positive = result_Ele_loss;
    result_Ele_loss_positive(result_Ele_loss_positive<0) = 0;
    result_Ele_loss_negtive = result_Ele_loss;
    result_Ele_loss_negtive(result_Ele_loss_negtive>0) = 0;
    
    stackedbar = @(x, A) bar(x, A, 'stacked');
    prettyline = @(x, y) plot(x, y, 'Color',firebrick, 'LineStyle','-','LineWidth',1.5, 'Marker', '.', 'MarkerSize', 13);
    
    figure
    if caseType == 31
        st = 2; en = 3;
    elseif caseType == 32
        %TODO
    else
        st = 1; en =3;
    end
    fig = 1;
    for IES_no = st : en
        subplot(en - st + 1, 1, fig)
        hold on
        bar_positive = [result_Ele_loss_positive(:,IES_no), result_CHP_power(:,IES_no), ...
            result_EH_solarP(:, IES_no) + result_EH_windP(:, IES_no), result_ES_discharge(:,IES_no)] ./ 1000;
        bar_negtive = [result_Ele_loss_negtive(:,IES_no), - result_eBoiler_E(:, IES_no), -result_ES_charge(:,IES_no)] ./1000;
        Egen = (result_Ele_loss(:,IES_no) + result_CHP_power(:,IES_no) + result_EH_solarP(:, IES_no) + result_EH_windP(:, IES_no)...
            + result_ES_discharge(:,IES_no) - result_ES_charge(:,IES_no)) ./1000 ;
        Eload = (result_EH_Le(:, IES_no) + result_EH_Edr(:, IES_no)) ./1000;
        Eload_base = result_EH_Le(:, IES_no) ./1000;
        
        yyaxis left;
        H1 = stackedbar(t1, bar_positive);
        H1(1).FaceColor = dodgerblue;
        H1(1).EdgeColor = 'none';
        H1(2).FaceColor = yellowgreen;
        H1(2).EdgeColor = 'none';
        H1(3).FaceColor = gold;
        H1(3).EdgeColor = 'none';
        H1(4).FaceColor = indianred;
        H1(4).EdgeColor = 'none';
        
        H3 = bar(bar_negtive, 'stacked');
        H3(1).FaceColor = H1(1).FaceColor;
        H3(1).EdgeColor = 'none';
        H3(3).FaceColor = H1(4).FaceColor;
        H3(3).EdgeColor = 'none';
        H3(2).FaceColor = chocolate3;
        H3(2).EdgeColor = 'none';
        H4 = plot(t1,[Eload_base, Eload]);
        set(H4(1), 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 13);
        set(H4(2), 'Color', 'black', 'LineStyle', '-', 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 13);

        ylabel('electric power(MW)');
        ylim([-1,2.5]);
        yticks(-1:0.5:2.5);
        
        yyaxis right;
        
        H2 = prettyline(t2, result_ES_SOC( : , IES_no)); 
        ylabel('SOC');
        ylim([0,1]);
        yticks(0.1:0.2:1);
        
        xlabel(sprintf('MES%d', IES_no));
        xlim([0, 24 * period + 1]);
        xticks(0:(24 * period / 4) : 24 * period);
        xticklabels({ '0:00','6:00','12:00','18:00','24:00' });
      
        if IES_no == st
            le = legend([H1(1), H1(2), H1(3), H1(4), H3(2), H2, H4(1), H4(2)],...
                'transformer','CHP','renewable energies','EES','EB','SOC of EES','base electric load','total electric load',...
                'Location','northoutside','Orientation','horizontal');
            set(le, 'Box', 'off');
        end
        set(gcf,'Position',[0 0 550 500]);
        fig = fig + 1;
    end
    %%%%%%%%%%%%%%
    
    figure
    if caseType == 31
        st = 3; en = 3;
    elseif caseType == 32
        %TODO
    else
        st = 1; en =3;
    end
    fig = 1;
    for IES_no = st : en
        subplot(en - st + 1, 1, fig)
        hold on
        bar_positive = [result_CHP_heat(:,IES_no), result_Boiler_heat(:,IES_no), result_eBoiler_H(:,IES_no), result_HS_discharge(:,IES_no)] ./1000;
        bar_negtive = -result_HS_charge(:,IES_no) ./1000;
        Hgen = (result_CHP_heat(:,IES_no) + result_Boiler_heat(:,IES_no) + result_HS_discharge(:,IES_no)...
            - result_HS_charge(:,IES_no) + result_eBoiler_H(:,IES_no)) ./1000;
        Hload = (result_EH_Lh(:, IES_no) + result_EH_Hdr(:, IES_no)) ./1000;
        Hload_base = result_EH_Lh(:, IES_no) ./1000;
        
        yyaxis left;
        H1 = stackedbar(t1, bar_positive);
        H1(1).FaceColor = yellowgreen;
        H1(1).EdgeColor = 'none';
        H1(2).FaceColor = gold;
        H1(2).EdgeColor = 'none';
        H1(3).FaceColor = chocolate3;
        H1(3).EdgeColor = 'none';
        H1(4).FaceColor = indianred;
        H1(4).EdgeColor = 'none';
        
        H3 = bar(bar_negtive,'stacked');
        H3(1).FaceColor = H1(4).FaceColor;
        H3(1).EdgeColor = 'none';
        
        H4 = plot(t1,[Hload_base, Hload]);
        set(H4(1), 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 13);
        set(H4(2), 'Color', 'black', 'LineStyle', '-', 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 13);
       
        ylabel('thermal power(MW)');
        ylim([-1, 2]);
        yticks(-1 : 0.5 : 2);
        
        yyaxis right;
        H2 = prettyline(t2, result_HS_SOC(:,IES_no));
        ylabel('SOC');
        ylim([0,1]);
        yticks(0.1:0.2:1);
        
        xlabel(sprintf('MES%d',IES_no));
        xlim([0, 24 * period + 1]);
        xticks(0 : (24 * period / 4) : 24 * period);
        xticklabels({'0:00','6:00','12:00','18:00','24:00'});

        if IES_no == st
            le = legend([H1(1),H1(2),H1(3),H1(4),H2,H4(1),H4(2)],...
                'CHP','GF','EB','TES','SOC of TES','base thermal load','total thermal load',...
                'Location','northoutside','Orientation','horizontal');
            set(le,'Box','off');
        end
        set(gcf,'Position',[0 0 500 500]);
        fig = fig + 1;
    end
end

if caseType == 32
minMarketPrice = 0.1;
maxMarketPrice = 1;
figure
    hold on;
     price = minMarketPrice : 0.01 : 1;
    
    H1= plot(price, demand(1, : ) / 1000);
    H3= plot(price, demandNoThss(1, : ) / 1000);
    H4= plot(price, demandNoSS(1, : ) / 1000);
    H2= plot(price, demandNoESS(1, : ) / 1000);
    
    set(H1(1), 'Color', 'red', 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 12);
    set(H2(1), 'Color', 'black', 'LineWidth',1.5);
    set(H3(1), 'Color', 'black', 'LineStyle', ':', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12);
    set(H4(1), 'Color', royalblue, 'LineStyle', ':', 'LineWidth', 3);
%     set(H4(2),'Color',yellowgreen, 'LineStyle',':','LineWidth',1.5,'Marker', '*', 'MarkerSize', 5);
%     set(H4(3),'Color',indianred, 'LineStyle',':','LineWidth',1.5,'Marker', '*', 'MarkerSize', 5);
%     set(H1(2),'Color',yellowgreen,'LineWidth',1.5, 'Marker', '.', 'MarkerSize', 13);
%     set(H1(3),'Color',indianred, 'LineWidth',1.5, 'Marker', '.', 'MarkerSize', 13);
%     set(H2(2),'Color',yellowgreen, 'LineStyle','-.','LineWidth',1.5);
%     set(H2(3),'Color',indianred, 'LineStyle','-.','LineWidth',1.5);
%     set(H3(2),'Color',yellowgreen,'LineWidth',1.5);
%     set(H3(3),'Color',indianred,'LineWidth',1.5);
    xlabel('electricity price（yuan)')
    ylabel('demand for electricity(MW)')
    le = legend([H1(1),H2(1),H3(1),H4(1)],...
        'EES+TES', 'TES','EES','no storages',...
        'Location','northoutside','Orientation','horizontal');
    set(le,'Box','off');
    set(gcf,'Position',[0 0 500 200]);

end

