clc; clear; close all
% 数据处理
load('../collaborate.mat');
result_Ele_collaborate = result_Ele; priceArray_collaborate = priceArray;
load('../autonomous.mat', 'result_Ele', 'priceArray');
result_Ele_autonomous = result_Ele; priceArray_autonomous = priceArray;
load('../collaborate_feedin.mat', 'result_Ele', 'priceArray');
result_Ele_collaborate_feedin = result_Ele; priceArray_collaborate_feedin = priceArray;
result_Ele = result_Ele_collaborate; priceArray = priceArray_collaborate;
global period
plotAux;
%----------------负荷和可再生能源的曲线----------------
    figure
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
        eval(['EH_Le_origin = EH_Le_base + EH', num2str(IES_no),'_Le_drP_total / sum(EH', num2str(IES_no),'_Le_flag) * EH', num2str(IES_no),'_Le_flag;']);
        eval(['EH_Lh_origin = EH_Lh_base + EH', num2str(IES_no),'_Lh_drP_total / sum(EH', num2str(IES_no),'_Lh_flag) * EH', num2str(IES_no),'_Lh_flag;']);

        subplot(3 , 2 , (IES_no - 1) * 2 + 1 )
        hold on;
        plot(t, EH_Le_base / 1000, 'Color', 'b', 'LineStyle', '-', 'LineWidth', w)
        plot(t, EH_Le_origin / 1000, 'Color', 'b', 'LineStyle', '-.', 'LineWidth', w) 
        plot(t, EH_Lh_base / 1000, 'Color', 'r', 'LineStyle', '-', 'LineWidth',w)
        plot(t, EH_Lh_origin / 1000, 'Color', 'r', 'LineStyle', '-.', 'LineWidth', w) 
        xlim([0, 24 * period])
        ylim([0, max(max(EH_Le), max(EH_Lh)) / 1000]);
        xticks(0 : (24 * period / 4) : 24 * period);
        xticklabels({ '0:00', '6:00', '12:00', '18:00', '24:00' });
        xlabel(sprintf('MES_%d', IES_no))
        ylabel('Load(MW)')
        if IES_no == 1
            H1 = legend('fixed electric load','total electric load','fixed thermal load','total thermal load',...
                'Location','northoutside','Orientation','horizontal');
            set(H1,'Box','off');
        end
 
        subplot(3,2,(IES_no - 1) * 2 + 2)
        hold on
        if solar_max(IES_no) > 0
            plot(t, EH_solarP / 1000, 'Color', 'k', 'LineStyle', '-', 'LineWidth', w);
            le = legend('PV','Location','northoutside','Orientation','horizontal');
        end
        if wind_max(IES_no) > 0
            plot(t, EH_windP / 1000, 'Color', 'k', 'LineStyle', '--', 'LineWidth', w);
            le = legend('wind','Location','northoutside','Orientation','horizontal');
        end
        set(le,'Box','off');
        xlim([0, 24*period]);
        xticks(0 : (24 * period /4) : 24 * period);
        xticklabels({'0:00','6:00','12:00','18:00','24:00'});
        xlabel(sprintf('MES_%d', IES_no))
        ylabel('power(MW)')
            
    end
    set(gcf,'Position',[0 0 650 500]);
    
    
    %------------------数据处理----------------
    result_Gas = result_CHP_G + result_Boiler_G;
    for IES_no = 1 : 3
        eval(['result_Ele_loss(:,IES_no) = result_Ele(:,IES_no) .* EH',num2str(IES_no),'.Ele_eff;']); % eleLimit(3)是线损率
        eval(['result_CHP_power(:,IES_no) = result_CHP_G(:,IES_no) .* EH',num2str(IES_no),'.CHP_GE_eff; ']);
        eval(['result_CHP_heat(:,IES_no) = result_CHP_G(:,IES_no) .* EH',num2str(IES_no),'.CHP_GH_eff; ']);
        eval(['result_Boiler_heat(:,IES_no) = result_Boiler_G(:,IES_no) .* EH',num2str(IES_no),'.Boiler_eff;']);
        eval(['result_eBoiler_H(:, IES_no) = result_eBoiler_E(:, IES_no) .* EH',num2str(IES_no),'.eBoiler_eff;']);
    end
    %------------------测试优化结果-------------
    ee = 1e-3;
    
    % 当储能的最大充、放电功率很大时，1000 * 1e-3 也会越线，因此应该提前将1e-3置为0
    result_ES_discharge(result_ES_discharge < ee) = 0;
    result_ES_charge(result_ES_charge < ee) = 0;
    result_HS_discharge(result_HS_discharge < ee) = 0;
    result_HS_charge(result_HS_charge < ee) = 0;
    
    % 计算
    % 计算总成本 按网价计算
    if exist('priceArray', 'var')
        cost_clear =  (result_Ele' * priceArray + sum(result_Gas)' * gasPrice1) / period;
    end
    totalCost = zeros(1,3);
    totalCost(1) = sum((result_Ele_autonomous' *  elePrice + sum(result_Gas)' * gasPrice1) / period);
    totalCost(2) = sum((result_Ele_collaborate_feedin' *  elePrice + sum(result_Gas)' * gasPrice1) / period);
    totalCost(3) = sum((result_Ele_collaborate' *  elePrice + sum(result_Gas)' * gasPrice1) / period);
    % ------------------绘图-------------------
    t1 = 1 : 1 : 24 * period;
    t2 = 0 : 1 : 24 * period;
    optNumber = 24;
    w=1.5;
    %----------------协同与不协同各MES对比----
    [c4_gridClearDemand_autonomous] = drawMES_stacked(t1, result_Ele_autonomous, EH_res_total, -3e3, 4e3, '../autonomous.mat', elePrice);
    [c4_gridClearDemand_collaborate_feedin] = drawMES_stacked(t1, result_Ele_collaborate_feedin, EH_res_total, -3e3, 4e3, '../collaborate_feedin.mat', elePrice);
    [c4_gridClearDemand_collaborate] = drawMES_stacked(t1, result_Ele_collaborate, EH_res_total, -3e3, 4e3, '../collaborate.mat', elePrice);

    %-----------------阻塞管理-----------------
    figure;
    hold on;
    H1 = plot(t1, [c4_gridClearDemand_autonomous, c4_gridClearDemand_collaborate_feedin, c4_gridClearDemand_collaborate, ]/1000);
    set(H1(1), 'Color', 'Black', 'LineWidth', 1.5, 'LineStyle', '--', 'Marker', '.', 'MarkerSize', 13);
    set(H1(2), 'Color', 'Black', 'LineWidth', 1.5);
    set(H1(3), 'Color', firebrick, 'LineWidth', 1.5);
    stairs(t2, ones(24*period+1, 1) .* 0, 'Color',gray,'LineStyle','--','LineWidth',1);
    stairs(t2, ones(24*period+1, 1) .* eleLimit_total(1)/1000, 'Color',gray,'LineStyle','--','LineWidth',1);
    stairs(t2, - ones(24*period+1, 1) * eleLimit_total(1)/1000, 'Color',gray,'LineStyle','--','LineWidth',1);
    ylabel('transformer power(MW)');
    uplimit= max(c4_gridClearDemand_collaborate + EH_res_total) / 1000 * 1.1;
    lowerlimit=-eleLimit_total(2) / 1000 * 1.1;
    
    le = legend([H1(1), H1(2), H1(3)],...
        'autonomous','collaborative autonomous(feed in)', 'collaborative autonomous',...
        'Orientation', 'vertical');...
    
    set(le,'Box','off');
    xlim([0, 25]); xticks(0 : 3: 24); xticklabels({'0:00','3:00', '6:00','9:00','12:00','15:00','18:00','21:00','24:00'});
    set(gcf,'Position',[0 0 500 200]);
    
    %--------------消纳率计算---------------
   accomodation = zeros(1,3);
   accomodation(1) = calculate_accomodation_rate(result_Ele_autonomous, c4_gridClearDemand_autonomous, ...
       EH_res_total, res_total, result_EH_windP, result_EH_solarP );
   accomodation(2) = calculate_accomodation_rate(result_Ele_collaborate_feedin, c4_gridClearDemand_collaborate_feedin,...
       EH_res_total, res_total, result_EH_windP, result_EH_solarP );
   accomodation(3) = calculate_accomodation_rate(result_Ele_collaborate, c4_gridClearDemand_collaborate,...
       EH_res_total, res_total, result_EH_windP, result_EH_solarP );
    totalCost
    accomodation
   %----------------优化结果2---------------
    result_Ele_loss_positive = result_Ele_loss;
    result_Ele_loss_positive(result_Ele_loss_positive<0) = 0;
    result_Ele_loss_negtive = result_Ele_loss;
    result_Ele_loss_negtive(result_Ele_loss_negtive>0) = 0;

    figure
    st = 1; en =3;
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
        H3(2).FaceColor = gray;
        H3(2).EdgeColor = 'none';
        H4 = plot(t1,[Eload_base, Eload]);
        set(H4(1), 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 1.5);
        set(H4(2), 'Color', 'black', 'LineStyle', '-', 'LineWidth', 1.5);

        ylabel('power(MW)');
        ylim([min(sum(bar_negtive, 2)) * 1.1 - 0.01 ,max(sum(bar_positive, 2)) * 1.1]);
        
        if max(result_ES_discharge(:, IES_no)) > 0 && max(result_ES_charge(:, IES_no)) > 0
            yyaxis right;
            H2 = prettyline(t2, result_ES_SOC( : , IES_no)); 
            ylabel('SOC');
            ylim([0,1]);
            yticks(0.1:0.2:1);
        end
        
        xlabel(sprintf('MES_%d', IES_no));
        xlim([0, 24 * period + 1]);
        xticks(0:(24 * period / 4) : 24 * period);
        xticklabels({ '0:00','6:00','12:00','18:00','24:00' });
      
        if max(result_ES_discharge(:, IES_no)) > 0 && max(result_ES_charge(:, IES_no)) > 0
            le = legend([H1(1), H1(2), H1(3), H1(4), H3(2), H2, H4(1), H4(2)],...
                'imported/exported electricity','CHP','renewable energies','EES','EB',...
                'SOC of EES','fixed electric load','total electric load',...
                'Location','northoutside','Orientation','horizontal');
            set(le, 'Box', 'off');
        end
%         set(gcf,'Position',[0 0 660 500]);
        set(le, 'NumColumns', 5);
        set(gcf,'Position',[0 0 590 500]);
        fig = fig + 1;
    end
    
    figure
    st = 1; en =3;
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
        H1(3).FaceColor = gray;
        H1(3).EdgeColor = 'none';
        H1(4).FaceColor = indianred;
        H1(4).EdgeColor = 'none';
        
        H3 = bar(bar_negtive,'stacked');
        H3(1).FaceColor = H1(4).FaceColor;
        H3(1).EdgeColor = 'none';
        
        H4 = plot(t1,[Hload_base, Hload]);
        set(H4(1), 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 1.5);
        set(H4(2), 'Color', 'black', 'LineStyle', '-', 'LineWidth', 1.5);
       
        ylabel('power(MW)');
        ylim([min(sum(bar_negtive, 2)) * 1.1 - 0.01,max(sum(bar_positive, 2)) * 1.1]);
        if max(result_HS_discharge(:, IES_no)) > 10 && max(result_HS_charge(:, IES_no)) > 10
            yyaxis right;
            H2 = prettyline(t2, result_HS_SOC(:,IES_no));
            ylabel('SOC');
            ylim([0,1]);
            yticks(0.1:0.2:1);
        end
       
        
        xlabel(sprintf('MES_%d',IES_no));
        xlim([0, 24 * period + 1]);
        xticks(0 : (24 * period / 4) : 24 * period);
        xticklabels({'0:00','6:00','12:00','18:00','24:00'});

        if max(result_HS_discharge(:, IES_no)) > 10 && max(result_HS_charge(:, IES_no)) > 10
            le = legend([H1(1),H1(2),H1(3),H1(4),H2,H4(1),H4(2)],...
                'CHP','GF','EB','TES',...
                'SOC of TES','fixed thermal load','total thermal load',...
                'Location','northoutside','Orientation','horizontal');
            set(le,'Box','off');
            set(le, 'NumColumns', 4);
        end
%         set(gcf,'Position',[0 0 660 500]);
        set(gcf,'Position',[0 0 590 500]);
        fig = fig + 1;
    end