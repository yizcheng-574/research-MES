% close all
% 数据处理
global period

%--------------------------------------负荷和可再生能源的曲线--------------------------------------

figure(1)
optNumber=24;
t=1:1:24*period;
w=1.2;

% EH1_Le(24*period+1) = EH1_Le(24*period);
% EH1_Lh(24*period+1) = EH1_Lh(24*period);
% EH1_solarP(24*period+1) = EH1_solarP(24*period);
% EH1_windP(24*period+1) = EH1_windP(24*period);
% EH2_Le(24*period+1) = EH2_Le(24*period);
% EH2_Lh(24*period+1) = EH2_Lh(24*period);
% EH2_solarP(24*period+1) = EH2_solarP(24*period);
% EH2_windP(24*period+1) = EH2_windP(24*period);
% EH3_Le(24*period+1) = EH3_Le(24*period);
% EH3_Lh(24*period+1) = EH3_Lh(24*period);
% EH3_solarP(24*period+1) = EH3_solarP(24*period);
% EH3_windP(24*period+1) = EH3_windP(24*period);


for IES_no = 1 : 3
    eval(['EH_Le = EH',num2str(IES_no),'_Le;']);
    eval(['EH_Lh = EH',num2str(IES_no),'_Lh;']);
    eval(['EH_solarP = EH',num2str(IES_no),'_solarP;']);
    eval(['EH_windP = EH',num2str(IES_no),'_windP;']);
    subplot(3 , 2 , (IES_no - 1) * 2 + 1 )
    hold on;
    stairs(t,EH_Le/1000,'Color','b','LineStyle','--','LineWidth',w) %plot
    stairs(t,EH_Lh/1000,'Color','r','LineStyle','-','LineWidth',w)
    xlim([0,24*period])
    set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',{'0:00','6:00','12:00','18:00','24:00'})
    ylabel('负荷需求(MW)')
    legend('电负荷','热负荷','Location','northoutside','Orientation','horizontal')
    % ylabel('load / kW')
    % xlabel('time / h')
    % legend('Le','Lh','Location','northoutside','Orientation','horizontal')
    % xlabel('时间(h)')
    subplot(3,2,(IES_no - 1) * 2 + 2 )
    hold on
    stairs(t,EH_solarP/1000,'Color','r','LineStyle','-','LineWidth',w)
    stairs(t,EH_windP/1000,'Color','b','LineStyle','--','LineWidth',w)
    xlim([0,24*period])
    set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',{'0:00','6:00','12:00','18:00','24:00'})
    % ylabel('RES power / kW')
    % xlabel('time / h')
    % legend('PV','WT','Location','northoutside','Orientation','horizontal')
    ylabel('发电功率(MW)')
    % xlabel('时间(h)')
    legend('光伏','风电','Location','northoutside','Orientation','horizontal')
end

%--------------------------------------数据处理--------------------------------------
result_Gas = result_CHP_G + result_Boiler_G;
for IES_no = 1 : 3
    eval(['result_Ele_loss(:,IES_no) = result_Ele(:,IES_no) .* eleLimit',num2str(IES_no),'(3);']); % eleLimit(3)是线损率
    eval(['result_CHP_power(:,IES_no) = result_CHP_G(:,IES_no) .* CHP',num2str(IES_no),'_para(1); ']);
    eval(['result_CHP_heat(:,IES_no) = result_CHP_G(:,IES_no) .* CHP',num2str(IES_no),'_para(2); ']); 
    eval(['result_Boiler_heat(:,IES_no) = result_Boiler_G(:,IES_no) .* Boiler',num2str(IES_no),'_para(1);']);  
end


%
%--------------------------------------测试优化结果--------------------------------------
% eleLimit(3)是线损率
ee = 1e-3;

% 当储能的最大充、放电功率很大时，1000 * 1e-3 也会越线，因此应该提前将1e-3置为0
result_ES_discharge(result_ES_discharge<ee) = 0;
result_ES_charge(result_ES_charge<ee) = 0;
result_HS_discharge(result_HS_discharge<ee) = 0;
result_HS_charge(result_HS_charge<ee) = 0;

% IES1
%电、热功率平衡性测试，降低要求至大于零
result_balance_P_1 = result_Ele_loss(:,1) + result_CHP_power(:,1) + result_ES_discharge(:,1) - result_ES_charge(:,1) - EH1_Le + EH1_windP + EH1_solarP;
result_balance_H_1 = result_CHP_heat(:,1) + result_Boiler_heat(:,1) + result_HS_discharge(:,1) - result_HS_charge(:,1) - EH1_Lh;
%充、放功率至少有一个是零
result_check_ES_1 = result_ES_discharge(:,1) .* result_ES_charge(:,1);
result_check_HS_1 = result_HS_discharge(:,1) .* result_HS_charge(:,1); %有一点小问题，因为热过于充裕
% 20180129 对充放电乘积约束的再思考
if max(result_check_ES_1) > 0
    disp(['EH1 电储能需要修正'])
    for i=1:length(result_check_ES_1)
        if result_check_ES_1(i) > 0
            if result_balance_P_1(i) < ee % 如果充放同时进行，而又没有多余的电，才算有问题
                disp(['EH1 电储能结果有问题 !!! 时间是', num2str(i)])
            else % 对充放同时的场景进行修正：SOC不变，充放电功率保留一个非零，且值有变化
                deltaSOC = result_ES_SOC(i+1,1) - result_ES_SOC(i,1);
                if deltaSOC > 0 % 表示充电
                    result_ES_charge(i,1) = deltaSOC * ES1_para(1) / ES1_para(7); % 乘以容量，除以效率
                    result_ES_discharge(i,1) = 0;
                else % 表示放电
                    result_ES_charge(i,1) = 0;
                    result_ES_discharge(i,1) = - deltaSOC * ES1_para(1) * ES1_para(7); % 乘以容量，乘以效率;
                end
            end
        end
    end
    % 完成修正后重新计算
    result_balance_P_1 = result_Ele_loss(:,1) + result_CHP_power(:,1) + result_ES_discharge(:,1) - result_ES_charge(:,1) - EH1_Le + EH1_windP + EH1_solarP;
    result_check_ES_1 = result_ES_discharge(:,1) .* result_ES_charge(:,1);
end
if max(result_check_HS_1) > 0
    disp(['EH1 热储能需要修正'])
    for i=1:length(result_check_HS_1)
        if result_check_HS_1(i) > 0
            if result_balance_H_1(i) < ee % 如果充放同时进行，而又没有多余的电，才算有问题
                disp(['EH1 热储能结果有问题 !!! 时间是', num2str(i)])
            else % 对充放同时的场景进行修正：SOC不变，充放电功率保留一个非零，且值有变化
                deltaSOC = result_HS_SOC(i+1,1) - result_HS_SOC(i,1);
                if deltaSOC > 0 % 表示充电
                    result_HS_charge(i,1) = deltaSOC * HS1_para(1) / HS1_para(7); % 乘以容量，除以效率
                    result_HS_discharge(i,1) = 0;
                else % 表示放电
                    result_HS_charge(i,1) = 0;
                    result_HS_discharge(i,1) = - deltaSOC * HS1_para(1) * HS1_para(7); % 乘以容量，乘以效率;
                end
            end
        end
    end
    % 完成修正后重新计算
    result_balance_H_1 = result_CHP_heat(:,1) + result_Boiler_heat(:,1) + result_HS_discharge(:,1) - result_HS_charge(:,1) - EH1_Lh;
    result_check_HS_1 = result_HS_discharge(:,1) .* result_HS_charge(:,1);
end



% IES2
%电、热功率平衡性测试，降低要求至大于零
result_balance_P_2 = result_Ele_loss(:,2) + result_CHP_power(:,2) + result_ES_discharge(:,2) - result_ES_charge(:,2) - EH2_Le + EH2_windP + EH2_solarP;
result_balance_H_2 = result_CHP_heat(:,2) + result_Boiler_heat(:,2) + result_HS_discharge(:,2) - result_HS_charge(:,2) - EH2_Lh;
%充、放功率至少有一个是零
result_check_ES_2 = result_ES_discharge(:,2) .* result_ES_charge(:,2);
result_check_HS_2 = result_HS_discharge(:,2) .* result_HS_charge(:,2); %有一点小问题，因为热过于充裕
% 20180129 对充放电乘积约束的再思考
if max(result_check_ES_2) > 0
    disp(['EH2 电储能需要修正'])
    for i=1:length(result_check_ES_2)
        if result_check_ES_2(i) > 0
            if result_balance_P_2(i) < ee % 如果充放同时进行，而又没有多余的电，才算有问题
                disp(['EH2 电储能结果有问题 !!! 时间是', num2str(i)])
            else % 对充放同时的场景进行修正：SOC不变，充放电功率保留一个非零，且值有变化
                deltaSOC = result_ES_SOC(i+1,2) - result_ES_SOC(i,2);
                if deltaSOC > 0 % 表示充电
                    result_ES_charge(i,2) = deltaSOC * ES2_para(1) / ES2_para(7); % 乘以容量，除以效率
                    result_ES_discharge(i,2) = 0;
                else % 表示放电
                    result_ES_charge(i,2) = 0;
                    result_ES_discharge(i,2) = - deltaSOC * ES2_para(1) * ES2_para(7); % 乘以容量，乘以效率;
                end
            end
        end
    end
    % 完成修正后重新计算
    result_balance_P_2 = result_Ele_loss(:,2) + result_CHP_power(:,2) + result_ES_discharge(:,2) - result_ES_charge(:,2) - EH2_Le + EH2_windP + EH2_solarP;
    result_check_ES_2 = result_ES_discharge(:,2) .* result_ES_charge(:,2);
end
if max(result_check_HS_2) > 0
    disp(['EH2 热储能需要修正'])
    for i=1:length(result_check_HS_2)
        if result_check_HS_2(i) > 0
            if result_balance_H_2(i) < ee % 如果充放同时进行，而又没有多余的电，才算有问题
                disp(['EH2 热储能结果有问题 !!! 时间是', num2str(i)])
            else % 对充放同时的场景进行修正：SOC不变，充放电功率保留一个非零，且值有变化
                deltaSOC = result_HS_SOC(i+1,2) - result_HS_SOC(i,2);
                if deltaSOC > 0 % 表示充电
                    result_HS_charge(i,2) = deltaSOC * HS2_para(1) / HS2_para(7); % 乘以容量，除以效率
                    result_HS_discharge(i,2) = 0;
                else % 表示放电
                    result_HS_charge(i,2) = 0;
                    result_HS_discharge(i,2) = - deltaSOC * HS2_para(1) * HS2_para(7); % 乘以容量，乘以效率;
                end
            end
        end
    end
    % 完成修正后重新计算
    result_balance_H_2 = result_CHP_heat(:,2) + result_Boiler_heat(:,2) + result_HS_discharge(:,2) - result_HS_charge(:,2) - EH2_Lh;
    result_check_HS_2 = result_HS_discharge(:,2) .* result_HS_charge(:,2);
end



% IES3
%电、热功率平衡性测试，降低要求至大于零
result_balance_P_3 = result_Ele_loss(:,3) + result_CHP_power(:,3) + result_ES_discharge(:,3) - result_ES_charge(:,3) - EH3_Le + EH3_windP + EH3_solarP;
result_balance_H_3 = result_CHP_heat(:,3) + result_Boiler_heat(:,3) + result_HS_discharge(:,3) - result_HS_charge(:,3) - EH3_Lh;
%充、放功率至少有一个是零
result_check_ES_3 = result_ES_discharge(:,3) .* result_ES_charge(:,3);
result_check_HS_3 = result_HS_discharge(:,3) .* result_HS_charge(:,3); %有一点小问题，因为热过于充裕
% 20180129 
if max(result_check_ES_3) > 0
    disp(['EH3 电储能需要修正'])
    for i=1:length(result_check_ES_3)
        if result_check_ES_3(i) > 0
            if result_balance_P_3(i) < ee % 如果充放同时进行，而又没有多余的电，才算有问题
                disp(['EH3 电储能结果有问题 !!! 时间是', num2str(i)])
            else % 对充放同时的场景进行修正：SOC不变，充放电功率保留一个非零，且值有变化
                deltaSOC = result_ES_SOC(i+1,3) - result_ES_SOC(i,3);
                if deltaSOC > 0 % 表示充电
                    result_ES_charge(i,3) = deltaSOC * ES3_para(1) / ES3_para(7); % 乘以容量，除以效率
                    result_ES_discharge(i,3) = 0;
                else % 表示放电
                    result_ES_charge(i,3) = 0;
                    result_ES_discharge(i,3) = - deltaSOC * ES3_para(1) * ES3_para(7); % 乘以容量，乘以效率;
                end
            end
        end
    end
    % 完成修正后重新计算
    result_balance_P_3 = result_Ele_loss(:,3) + result_CHP_power(:,3) + result_ES_discharge(:,3) - result_ES_charge(:,3) - EH3_Le + EH3_windP + EH3_solarP;
    result_check_ES_3 = result_ES_discharge(:,3) .* result_ES_charge(:,3);
end
if max(result_check_HS_3) > 0
    disp(['EH3 热储能需要修正'])
    for i=1:length(result_check_HS_3)
        if result_check_HS_3(i) > 0
            if result_balance_H_3(i) < ee % 如果充放同时进行，而又没有多余的电，才算有问题
                disp(['EH3 热储能结果有问题 !!! 时间是', num2str(i)])
            else % 对充放同时的场景进行修正：SOC不变，充放电功率保留一个非零，且值有变化
                deltaSOC = result_HS_SOC(i+1,3) - result_HS_SOC(i,3);
                if deltaSOC > 0 % 表示充电
                    result_HS_charge(i,3) = deltaSOC * HS3_para(1) / HS3_para(7); % 乘以容量，除以效率
                    result_HS_discharge(i,3) = 0;
                else % 表示放电
                    result_HS_charge(i,3) = 0;
                    result_HS_discharge(i,3) = - deltaSOC * HS3_para(1) * HS3_para(7); % 乘以容量，乘以效率;
                end
            end
        end
    end
    % 完成修正后重新计算
    result_balance_H_3 = result_CHP_heat(:,3) + result_Boiler_heat(:,3) + result_HS_discharge(:,3) - result_HS_charge(:,3) - EH3_Lh;
    result_check_HS_3 = result_HS_discharge(:,3) .* result_HS_charge(:,3); %有一点小问题，因为热过于充裕
end



%检查馈线的平衡性，并检查馈线的功率是否越限
result_balance_grid = gridClearDemand + sum(result_Ele,2); %2表示按列相加

%计算总的购电功率和总的购气功率(查看是否越限)
% totalE = result_Ele(:,1) + result_Ele(:,2) + result_Ele(:,3);
% totalGas = EH1_G_CHP + EH1_G_Boiler + EH2_G_CHP + EH2_G_Boiler

%}
% 计算
%计算总成本 按网价计算
for IES_no = 1 : 3
    eval(['totalCost',num2str(IES_no),' = ( sum(result_Ele(:,',num2str(IES_no),') .* elePrice) + sum(result_Gas(:,',num2str(IES_no),') .* gasPrice',num2str(IES_no),') ) / period;']);
end
disp(['IES1总成本为 ',num2str(totalCost1),' 元'])
disp(['IES2总成本为 ',num2str(totalCost2),' 元'])
disp(['IES3总成本为 ',num2str(totalCost3),' 元'])
disp(['总成本为 ',num2str(totalCost1 + totalCost2 + totalCost3),' 元'])

%按出清价格 计算成本
% totalCost1 = ( sum(result_Ele(:,1) .* priceArray) + sum(result_Gas(:,1) .* gasPrice1) ) / period; %priceArray保留的是日内出清价格
% totalCost2 = ( sum(result_Ele(:,2) .* priceArray) + sum(result_Gas(:,2) .* gasPrice1) ) / period;
% totalCost3 = ( sum(result_Ele(:,3) .* priceArray) + sum(result_Gas(:,3) .* gasPrice3) ) / period;
% 
% disp(['IES1总成本2为 ',num2str(totalCost1),' 元'])
% disp(['IES2总成本2为 ',num2str(totalCost2),' 元'])
% disp(['IES3总成本2为 ',num2str(totalCost3),' 元'])
% disp(['总成本2为 ',num2str(totalCost1 + totalCost2 + totalCost3),' 元'])


%计算弃风光率
waste_power1 = sum(result_balance_P_1) / sum(EH1_solarP + EH1_windP) * 100; % 分子分母的period抵消了
waste_power2 = sum(result_balance_P_2) / sum(EH2_solarP + EH2_windP) * 100;
waste_power3 = sum(result_balance_P_3) / sum(EH3_solarP + EH3_windP) * 100;
disp(['IES1弃风光率为 ',num2str(waste_power1),' %'])
disp(['IES2弃风光率为 ',num2str(waste_power2),' %'])
disp(['IES3弃风光率为 ',num2str(waste_power3),' %'])

%计算热浪费率
waste_heat1 = sum(result_balance_H_1) / sum(EH1_Lh) * 100; % 分子分母的period抵消了
waste_heat2 = sum(result_balance_H_2) / sum(EH2_Lh) * 100;
waste_heat3 = sum(result_balance_H_3) / sum(EH3_Lh) * 100;
disp(['IES1热浪费率为 ',num2str(waste_heat1),' %'])
disp(['IES2热浪费率为 ',num2str(waste_heat2),' %'])
disp(['IES3热浪费率为 ',num2str(waste_heat3),' %'])

% %计算迭代次数
% times_dayahead = sum(iterationTimes(:,1)) / length(iterationTimes(:,1));
% times_inday = sum(iterationTimes(:,2)) / length(iterationTimes(:,2));
% disp(['日前优化平均次数 ',num2str(times_dayahead)])
% disp(['日内优化平均次数 ',num2str(times_inday)])





% --------------------------------------绘图--------------------------------------
t1 = 1:1:24*period;
t2 = 0:1:24*period;
optNumber = 24;
w=1.2;

% IES内用电 用气
% 电：CHP，电储能，联络线；
% 热：CHP，锅炉，热储能
%--------------------------------------数据拓展--------------------------------------
% result_Ele_loss(24*period+1, :) = result_Ele_loss(24*period, :);
% result_CHP_power(24*period+1, :) = result_CHP_power(24*period, :);
% result_CHP_heat(24*period+1, :) = result_CHP_heat(24*period, :);
% result_Boiler_heat(24*period+1, :) = result_Boiler_heat(24*period, :);

%{
figure(4)
subplot(2,1,1)
hold on
[AX,H1,H2] = plotyy(t2, result_Ele_loss(:,1), t2, result_ES_SOC(:,1),'stairs','plot');
stairs(t2, result_CHP_power(:,1), 'Color','r','LineStyle','--','LineWidth',w);
% legend('lateral','CHP','EESS','Location','northoutside','Orientation','horizontal')
legend('支线交换功率','CHP电功率','ESS SOC','Location','northoutside','Orientation','horizontal')
set(H1,'Color','b','LineStyle',':','LineWidth',w)
set(H2,'Color',[101, 147, 74]./255,'LineWidth',w)
% set(get(AX(1),'Ylabel'),'String','electricity power / kW') 
set(get(AX(1),'Ylabel'),'String','电功率(kW)') 
set(get(AX(2),'Ylabel'),'String','SOC') 
% xlabel('time / h') 
% text(1,1,'(c)支线2')
set(AX(1),'xlim',[0,24*period])
set(AX(2),'xlim',[0,24*period])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',0:(optNumber/4):optNumber)
set(AX(2),'XTick',[],'XTickLabel',[])
set(AX(1),'ylim',[-100,650])
set(AX(1),'YTick',-100:150:650)
set(AX(2),'ylim',[0,1])
set(AX(2),'YTick',0:0.2:1)

subplot(2,1,2)
hold on
[AX,H1,H2] = plotyy(t2, result_Boiler_heat(:,1), t2, result_HS_SOC(:,1),'stairs','plot');
stairs(t2, result_CHP_heat(:,1), 'Color','r','LineStyle','--','LineWidth',w)
legend('GF热功率','CHP热功率','ThSS SOC','Location','northoutside','Orientation','horizontal')
set(H1,'Color','b','LineStyle',':','LineWidth',w)
set(H2,'Color',[101, 147, 74]./255,'LineWidth',w)
% set(get(AX(1),'Ylabel'),'String','thermal power / kW') 
set(get(AX(1),'Ylabel'),'String','热功率(kW)') 
set(get(AX(2),'Ylabel'),'String','SOC') 
% xlabel('time / h') 
xlabel('时间(h)') 
% text(1,1,'(c)支线2')
set(AX(1),'xlim',[0,24*period])
set(AX(2),'xlim',[0,24*period])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',0:(optNumber/4):optNumber)
set(AX(2),'XTick',[],'XTickLabel',[])
set(AX(1),'ylim',[0,750])
set(AX(1),'YTick',0:150:750)
set(AX(2),'ylim',[0,1])
set(AX(2),'YTick',0:0.2:1)

set(gcf,'Position',[0 0 400 300]);

figure(5)

subplot(2,1,1)
hold on
[AX,H1,H2] = plotyy(t2, result_Ele_loss(:,2), t2, result_ES_SOC(:,2),'stairs','plot');
stairs(t2, result_CHP_power(:,2), 'Color','r','LineStyle','--','LineWidth',w);
% legend('lateral','CHP','EESS','Location','northoutside','Orientation','horizontal')
legend('支线交换功率','CHP电功率','ESS SOC','Location','northoutside','Orientation','horizontal')
set(H1,'Color','b','LineStyle',':','LineWidth',w)
set(H2,'Color',[101, 147, 74]./255,'LineWidth',w)
% set(get(AX(1),'Ylabel'),'String','electricity power / kW') 
set(get(AX(1),'Ylabel'),'String','电功率(kW)') 
set(get(AX(2),'Ylabel'),'String','SOC') 
% xlabel('时间') 
% text(1,1,'(c)支线2')
set(AX(1),'xlim',[0,24*period])
set(AX(2),'xlim',[0,24*period])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',0:(optNumber/4):optNumber)
set(AX(2),'XTick',[],'XTickLabel',[])
set(AX(1),'ylim',[-100,400])
set(AX(1),'YTick',-100:100:400)
set(AX(2),'ylim',[0,1])
set(AX(2),'YTick',0:0.2:1)

subplot(2,1,2)
hold on
[AX,H1,H2] = plotyy(t2, result_Boiler_heat(:,2), t2, result_HS_SOC(:,2),'stairs','plot');
stairs(t2, result_CHP_heat(:,2), 'Color','r','LineStyle','--','LineWidth',w)
legend('GF热功率','CHP热功率','ThSS SOC','Location','northoutside','Orientation','horizontal')
set(H1,'Color','b','LineStyle',':','LineWidth',w)
set(H2,'Color',[101, 147, 74]./255,'LineWidth',w)
% set(get(AX(1),'Ylabel'),'String','thermal power / kW')  
set(get(AX(1),'Ylabel'),'String','热功率(kW)') 
set(get(AX(2),'Ylabel'),'String','SOC') 
% xlabel('time / h') 
xlabel('时间(h)') 
% text(1,1,'(c)支线2')
set(AX(1),'xlim',[0,24*period])
set(AX(2),'xlim',[0,24*period])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',0:(optNumber/4):optNumber)
set(AX(2),'XTick',[],'XTickLabel',[])
set(AX(1),'ylim',[0,400])
set(AX(1),'YTick',0:80:400)
set(AX(2),'ylim',[0,1])
set(AX(2),'YTick',0:0.2:1)


set(gcf,'Position',[0 0 400 300]);



figure(6)

subplot(2,1,1)
hold on
[AX,H1,H2] = plotyy(t2, result_Ele_loss(:,3), t2, result_ES_SOC(:,3),'stairs','plot');
stairs(t2, result_CHP_power(:,3), 'Color','r','LineStyle','--','LineWidth',w);
% legend('lateral','CHP','EESS','Location','northoutside','Orientation','horizontal')
legend('支线交换功率','CHP电功率','ESS SOC','Location','northoutside','Orientation','horizontal')
set(H1,'Color','b','LineStyle',':','LineWidth',w)
set(H2,'Color',[101, 147, 74]./255,'LineWidth',w)
% set(get(AX(1),'Ylabel'),'String','electricity power / kW') 
set(get(AX(1),'Ylabel'),'String','电功率(kW)') 
set(get(AX(2),'Ylabel'),'String','SOC') 
% xlabel('时间') 
% text(1,1,'(c)支线2')
set(AX(1),'xlim',[0,24*period])
set(AX(2),'xlim',[0,24*period])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',0:(optNumber/4):optNumber)
set(AX(2),'XTick',[],'XTickLabel',[])
set(AX(1),'ylim',[-60,240])
set(AX(1),'YTick',-60:60:240)
set(AX(2),'ylim',[0,1])
set(AX(2),'YTick',0:0.2:1)

subplot(2,1,2)
hold on
[AX,H1,H2] = plotyy(t2, result_Boiler_heat(:,3), t2, result_HS_SOC(:,3),'stairs','plot');
stairs(t2, result_CHP_heat(:,3), 'Color','r','LineStyle','--','LineWidth',w)
legend('GF热功率','CHP热功率','ThSS SOC','Location','northoutside','Orientation','horizontal')
set(H1,'Color','b','LineStyle',':','LineWidth',w)
set(H2,'Color',[101, 147, 74]./255,'LineWidth',w)
% set(get(AX(1),'Ylabel'),'String','thermal power / kW') 
set(get(AX(1),'Ylabel'),'String','热功率(kW)') 
set(get(AX(2),'Ylabel'),'String','SOC') 
% xlabel('time / h') 
xlabel('时间(h)') 
% text(1,1,'(c)支线2')
set(AX(1),'xlim',[0,24*period])
set(AX(2),'xlim',[0,24*period])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',0:(optNumber/4):optNumber)
set(AX(2),'XTick',[],'XTickLabel',[])
set(AX(1),'ylim',[0,400])
set(AX(1),'YTick',0:80:400)
set(AX(2),'ylim',[0,1])
set(AX(2),'YTick',0:0.2:1)


set(gcf,'Position',[0 0 400 300]);
%}

% -------只记录正常运行的出清电价---------
%{
elePrice(24*period+1, :) = elePrice(24*period, :);
priceArray_record(24*period+1, :) = priceArray_record(24*period, :);

figure(8)
w=1.2;

hold on
stairs(t2, elePrice, 'Color','b','LineStyle','--','LineWidth',w);
stairs(t2, priceArray_record(:,1), 'Color','r','LineStyle','-','LineWidth',w);
stairs(t2, priceArray_record(:,2), 'Color','k','LineStyle',':','LineWidth',w);

% legend('utility price','day-ahead clearing price','real-time clearing price','real-time clearing price（IES3 isolated）')
legend('电网电价','日前出清电价','日内出清电价')
axis([0 24*period 0.2 1.2])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',0:(optNumber/4):optNumber)
% xlabel('time / h')
xlabel('时间(h)') 
% ylabel('price / yuan/kWh')
ylabel('电价(元/kWh)') 
% set(gca,'FontSize',14) % 设置文字大小，同时影响坐标轴标注、图例、标题等。
% set(get(gca,'XLabel'),'FontSize',14);
% set(get(gca,'YLabel'),'FontSize',14);
set(gcf,'Position',[0 0 400 200]);
set(gca,'Fontname','Monospaced');
%}


% -------同时记录正常运行和突然离网时的出清电价---------
%{
priceArray_normal = [0.370008142338676;0.385717522693708;0.363534521676111;0.368375205528558;0.353264238197299;0.353915728172134;0.437959326413856;0.392929250442548;1.03292467099294;0.561401875187474;0.537300862415980;0.544080498984673;0.535957976900646;0.546204794157670;0.540864819374419;0.544925848460731;0.555158632619513;0.542910592622054;0.553443592092666;0.561598149278862;0.713071708347918;0.585047037404080;0.593975061642661;0.578366550751511;0.578366550751511];
elePrice(25,:) = elePrice(24,:);
priceArray_record(25,:) = priceArray_record(24,:);
priceArray_normal(25,:) = priceArray_normal(24,:);

figure(7)
w=1.2;

hold on
stairs(t2, elePrice, 'Color','k','LineStyle',':','LineWidth',w);
stairs(t2, priceArray_record(:,1), 'Color','b','LineStyle','-.','LineWidth',w);
stairs(t2, priceArray_normal, 'Color',[101, 147, 74]./255,'LineStyle','-','LineWidth',w);
stairs(t2, priceArray_record(:,2), 'Color','r','LineStyle','--','LineWidth',w);

% legend('utility price','day-ahead clearing price','real-time clearing price','real-time clearing price（IES3 isolated）')
legend('电网电价','日前出清电价','日内出清电价','IES意外离网后的日内出清电价')
axis([0 24 0.2 1.5])
set(gca,'XTick',0:(optNumber/4):optNumber)
% xlabel('time / h')
xlabel('时间(h)') 
% ylabel('price / yuan/kWh')
ylabel('电价(元/kWh)') 
% set(gca,'FontSize',14) % 设置文字大小，同时影响坐标轴标注、图例、标题等。
% set(get(gca,'XLabel'),'FontSize',14);
% set(get(gca,'YLabel'),'FontSize',14);
set(gcf,'Position',[0 0 400 300]);
%}




close all
w=1.5;
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

%--------------------阻塞管理---------------------

c4_clearingPrice = priceArray;
c4_gridClearDemand = gridClearDemand;
c2_gridClearDemand = [-15421.1011934133;-15281.3776317716;-15173.9791241128;-15326.6973919484;-15738.6978179487;-14870.0142607562;-14244.5901637662;-14698.0091062568;4824.99999999979;4824.99999515603;4824.99999991426;4824.99997642087;2432.68540315238;3624.52057098259;3593.34530132161;3544.67468690340;3424.55714773556;4824.99999996666;4824.99999999012;4824.99999999137;4824.99999994970;1044.62888886440;1604.21144959550;2476.39814681124];
                    
                     
                     
% elePrice(24*period+1, :) = elePrice(24*period, :);
% c4_clearingPrice(24*period+1, :) = c4_clearingPrice(24*period, :);
% c4_gridClearDemand(24*period+1, :) = c4_gridClearDemand(24*period, :);
% c2_gridClearDemand(24*period+1, :) = c2_gridClearDemand(24*period, :);

figure(11)
hold on
[AX,H1,H2] = plotyy(t1, [-c2_gridClearDemand/1000 , -c4_gridClearDemand/1000], t1, (c4_clearingPrice-elePrice), 'bar', 'plot');
H1(1).EdgeColor = dodgerblue;
H1(1).FaceColor = dodgerblue;
H1(2).EdgeColor = yellowgreen;
H1(2).FaceColor = yellowgreen;
set(H2,'Color',firebrick, 'LineStyle','-','LineWidth',1.5, 'Marker', '.', 'MarkerSize', 13)

set(get(AX(1),'Ylabel'),'String','馈线购电功率(MW)') 
set(get(AX(2),'Ylabel'),'String','电价偏移量(元/kWh)') 
% xlabel('时间(h)') 
set(AX(1),'xlim',[0,24*period+1])
set(AX(2),'xlim',[0,24*period+1])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',{'0:00','6:00','12:00','18:00','24:00'})
set(AX(2),'XTick',[],'XTickLabel',[])
set(AX(1),'ylim',[-17, 17])
set(AX(1),'YTick',-16:8:16)
set(AX(2),'ylim',[-0.53125, 0.53125])
set(AX(2),'YTick',-0.5:0.25:0.5)

H3 = stairs(t2, ones(24*period+1, 1).*16, 'Color',gray,'LineStyle','--','LineWidth',1);
stairs(t2, ones(24*period+1, 1).*(-16/4), 'Color',gray,'LineStyle','--','LineWidth',1);

H10 = legend([H1(1),H1(2),H2,H3],'方案2馈线功率','方案6馈线功率','方案6电价偏移量','馈线功率约束');
set(H10,'Box','off');

set(gcf,'Position',[0 0 400 200]);




%--------------------优化结果2---------------------
result_Ele_loss_positive = result_Ele_loss;
result_Ele_loss_positive(result_Ele_loss_positive<0) = 0;
result_Ele_loss_negtive = result_Ele_loss;
result_Ele_loss_negtive(result_Ele_loss_negtive>0) = 0;

stackedbar = @(x, A) bar(x, A, 'stacked');
prettyline = @(x, y) plot(x, y, 'Color',firebrick, 'LineStyle','-','LineWidth',1.5, 'Marker', '.', 'MarkerSize', 13);


figure(12)
subplot(3,1,1)
hold on

bar_positive = [result_Ele_loss_positive(:,1), result_CHP_power(:,1), EH1_solarP+EH1_windP, result_ES_discharge(:,1)] ./ 1000;
bar_negtive = [result_Ele_loss_negtive(:,1), -result_ES_charge(:,1)] ./1000;

[AX,H1,H2] = plotyy(t1, bar_positive, t2, result_ES_SOC(:,1),stackedbar,prettyline);
H1(1).EdgeColor = dodgerblue;
H1(1).FaceColor = H1(1).EdgeColor;
H1(2).EdgeColor = yellowgreen;
H1(2).FaceColor = H1(2).EdgeColor;
H1(3).EdgeColor = gold;
H1(3).FaceColor = H1(3).EdgeColor;
H1(4).EdgeColor = indianred;
H1(4).FaceColor = H1(4).EdgeColor;

H3 = bar(bar_negtive,'stacked');
H3(1).EdgeColor = H1(1).EdgeColor;
H3(1).FaceColor = H3(1).EdgeColor;
H3(2).EdgeColor = H1(4).EdgeColor;
H3(2).FaceColor = H3(2).EdgeColor;

set(get(AX(1),'Ylabel'),'String','电功率(MW)') 
set(get(AX(2),'Ylabel'),'String','SOC')
xlabel('(a) IES1') 
set(AX(1),'xlim',[0,24*period+1])
set(AX(2),'xlim',[0,24*period+1])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',{'0:00','6:00','12:00','18:00','24:00'}) % 0:(optNumber/4):optNumber)
set(AX(2),'XTick',[],'XTickLabel',[])
set(AX(1),'ylim',[-5,35])
set(AX(1),'YTick',-5:10:35)
set(AX(2),'ylim',[0.1,0.9])
set(AX(2),'YTick',0.1:0.2:0.9)

H10 = legend([H1(1),H1(2),H1(3),H1(4),H2], '支线购电','CHP','可再生能源','ESS放电','SOC','Location','northoutside','Orientation','horizontal');
set(H10,'Box','off');

subplot(3,1,2)
hold on

bar_positive = [result_Ele_loss_positive(:,2), result_CHP_power(:,2), EH2_solarP+EH2_windP, result_ES_discharge(:,2)] ./1000;
bar_negtive = [result_Ele_loss_negtive(:,2), -result_ES_charge(:,2)] ./1000;

[AX,H1,H2] = plotyy(t1, bar_positive, t2, result_ES_SOC(:,2),stackedbar,prettyline);
H1(1).EdgeColor = dodgerblue;
H1(1).FaceColor = H1(1).EdgeColor;
H1(2).EdgeColor = yellowgreen;
H1(2).FaceColor = H1(2).EdgeColor;
H1(3).EdgeColor = gold;
H1(3).FaceColor = H1(3).EdgeColor;
H1(4).EdgeColor = indianred;
H1(4).FaceColor = H1(4).EdgeColor;

H3 = bar(bar_negtive,'stacked');
H3(1).EdgeColor = H1(1).EdgeColor;
H3(1).FaceColor = H3(1).EdgeColor;
H3(2).EdgeColor = H1(4).EdgeColor;
H3(2).FaceColor = H3(2).EdgeColor;

set(get(AX(1),'Ylabel'),'String','电功率(MW)') 
set(get(AX(2),'Ylabel'),'String','SOC') 
xlabel('(b) IES2') 
set(AX(1),'xlim',[0,24*period+1])
set(AX(2),'xlim',[0,24*period+1])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',{'0:00','6:00','12:00','18:00','24:00'})
set(AX(2),'XTick',[],'XTickLabel',[])
set(AX(1),'ylim',[-0.25,1.75])
set(AX(1),'YTick',-0.25:0.5:1.75)
set(AX(2),'ylim',[0.1,0.9])
set(AX(2),'YTick',0.1:0.2:0.9)


subplot(3,1,3)
hold on

bar_positive = [result_Ele_loss_positive(:,3), result_CHP_power(:,3), EH3_solarP+EH3_windP, result_ES_discharge(:,3)] ./1000;
bar_negtive = [result_Ele_loss_negtive(:,3), -result_ES_charge(:,3)] ./1000;

[AX,H1,H2] = plotyy(t1, bar_positive, t2, result_ES_SOC(:,3),stackedbar,prettyline);
H1(1).EdgeColor = dodgerblue;
H1(1).FaceColor = H1(1).EdgeColor;
H1(2).EdgeColor = yellowgreen;
H1(2).FaceColor = H1(2).EdgeColor;
H1(3).EdgeColor = gold;
H1(3).FaceColor = H1(3).EdgeColor;
H1(4).EdgeColor = indianred;
H1(4).FaceColor = H1(4).EdgeColor;

H3 = bar(bar_negtive,'stacked');
H3(1).EdgeColor = H1(1).EdgeColor;
H3(1).FaceColor = H3(1).EdgeColor;
H3(2).EdgeColor = H1(4).EdgeColor;
H3(2).FaceColor = H3(2).EdgeColor;

set(get(AX(1),'Ylabel'),'String','电功率(MW)') 
set(get(AX(2),'Ylabel'),'String','SOC') 
xlabel('(c) IES3') 
% xlabel({'时间(h)';'(c) IES3'}) 
set(AX(1),'xlim',[0,24*period+1])
set(AX(2),'xlim',[0,24*period+1])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',{'0:00','6:00','12:00','18:00','24:00'})
set(AX(2),'XTick',[],'XTickLabel',[])
set(AX(1),'ylim',[-1, 7])
set(AX(1),'YTick',-1:2:7)
set(AX(2),'ylim',[0.1,0.9])
set(AX(2),'YTick',0.1:0.2:0.9)

set(gcf,'Position',[0 0 400 500]);






figure(13)

subplot(3,1,1)
hold on
bar_positive = [result_CHP_heat(:,1), result_Boiler_heat(:,1), result_HS_discharge(:,1)] ./1000;
bar_negtive = [-result_HS_charge(:,1)] ./1000;
[AX,H1,H2] = plotyy(t1, bar_positive, t2, result_HS_SOC(:,1),stackedbar,prettyline);
H1(1).EdgeColor = yellowgreen;
H1(1).FaceColor = H1(1).EdgeColor;
H1(2).EdgeColor = gold;
H1(2).FaceColor = H1(2).EdgeColor;
H1(3).EdgeColor = indianred;
H1(3).FaceColor = H1(3).EdgeColor;

H3 = bar(bar_negtive,'stacked');
H3(1).EdgeColor = H1(3).EdgeColor;
H3(1).FaceColor = H3(1).EdgeColor;

set(get(AX(1),'Ylabel'),'String','热功率(MW)') 
set(get(AX(2),'Ylabel'),'String','SOC') 
xlabel('(a) IES1') 
set(AX(1),'xlim',[0,24*period+1])
set(AX(2),'xlim',[0,24*period+1])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',{'0:00','6:00','12:00','18:00','24:00'})
set(AX(2),'XTick',[],'XTickLabel',[])
set(AX(1),'ylim',[0,40])
set(AX(1),'YTick',0:10:40)
set(AX(2),'ylim',[0.1,0.9])
set(AX(2),'YTick',0.1:0.2:0.9)

H10 = legend([H1(1),H1(2),H1(3),H2], 'CHP','GF','ThSS放电','SOC','Location','northoutside','Orientation','horizontal');
set(H10,'Box','off');


subplot(3,1,2)
hold on
bar_positive = [result_CHP_heat(:,2), result_Boiler_heat(:,2), result_HS_discharge(:,2)] ./1000;
bar_negtive = [-result_HS_charge(:,2)] ./1000;
[AX,H1,H2] = plotyy(t1, bar_positive, t2, result_HS_SOC(:,2),stackedbar,prettyline);
H1(1).EdgeColor = yellowgreen;
H1(1).FaceColor = H1(1).EdgeColor;
H1(2).EdgeColor = gold;
H1(2).FaceColor = H1(2).EdgeColor;
H1(3).EdgeColor = indianred;
H1(3).FaceColor = H1(3).EdgeColor;

H3 = bar(bar_negtive,'stacked');
H3(1).EdgeColor = H1(3).EdgeColor;
H3(1).FaceColor = H3(1).EdgeColor;

set(get(AX(1),'Ylabel'),'String','热功率(MW)') 
set(get(AX(2),'Ylabel'),'String','SOC') 
xlabel('(b) IES2') 
set(AX(1),'xlim',[0,24*period+1])
set(AX(2),'xlim',[0,24*period+1])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',{'0:00','6:00','12:00','18:00','24:00'})
set(AX(2),'XTick',[],'XTickLabel',[])
set(AX(1),'ylim',[-0.2, 1.0])
set(AX(1),'YTick',-0.2:0.3:1.0)
set(AX(2),'ylim',[0.1,0.9])
set(AX(2),'YTick',0.1:0.2:0.9)

subplot(3,1,3)
hold on
bar_positive = [result_CHP_heat(:,3), result_Boiler_heat(:,3), result_HS_discharge(:,3)] ./1000;
bar_negtive = [-result_HS_charge(:,3)] ./1000;
[AX,H1,H2] = plotyy(t1, bar_positive, t2, result_HS_SOC(:,3),stackedbar,prettyline);
H1(1).EdgeColor = yellowgreen;
H1(1).FaceColor = H1(1).EdgeColor;
H1(2).EdgeColor = gold;
H1(2).FaceColor = H1(2).EdgeColor;
H1(3).EdgeColor = indianred;
H1(3).FaceColor = H1(3).EdgeColor;

H3 = bar(bar_negtive,'stacked');
H3(1).EdgeColor = H1(3).EdgeColor;
H3(1).FaceColor = H3(1).EdgeColor;

set(get(AX(1),'Ylabel'),'String','热功率(MW)') 
set(get(AX(2),'Ylabel'),'String','SOC') 
xlabel('(c) IES3') 
% xlabel('时间(h)') 
set(AX(1),'xlim',[0,24*period+1])
set(AX(2),'xlim',[0,24*period+1])
set(gca,'XTick',0:(24*period/4):24*period, 'XTickLabel',{'0:00','6:00','12:00','18:00','24:00'})
set(AX(2),'XTick',[],'XTickLabel',[])
set(AX(1),'ylim',[-1.5, 6.5])
set(AX(1),'YTick',-1.5:2:6.5)
set(AX(2),'ylim',[0.1,0.9])
set(AX(2),'YTick',0.1:0.2:0.9)

set(gcf,'Position',[0 0 400 500]);

