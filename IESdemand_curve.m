% ES_totalC_in, ES_maxSOC_in, ES_minSOC_in, ES_currentSOC_in, ES_targetSOC_in, ES_chargeTime, ES_eff_in
        ES1_para = [1, 0.8, 0.2, 0.4, 0.4, 6, 0.9];
        ES2_para = [1, 0.85, 0.15, 0.4, 0.4, 6, 0.9];
        ES3_para = [1, 0.85, 0.15, 0.4, 0.4, 6, 0.9];
        % HS_totalC_in, HS_maxSOC_in, HS_minSOC_in, HS_currentSOC_in, HS_targetSOC_in, HS_chargeTime, HS_eff_in
        HS1_para = [1, 0.9, 0.1, 0.5, 0.5, 5, 0.9];
        HS2_para = [1, 0.9, 0.1, 0.5, 0.5, 5, 0.9];
        HS3_para = [1, 0.9, 0.1, 0.5, 0.5, 5, 0.9];EH1 = EH_local_170828_v3(eleLimit1, gasLimit1, EH1_Le, EH1_Lh, EH1_solarP, EH1_windP, CHP1_para, Boiler1_para, ES1_para, HS1_para, dev_L, dev_PV, dev_WT, EH1_solarP_rate, EH1_windP_rate, EH1_Le_drP_rate, EH1_Le_drP_total, EH1_Lh_drP_rate, EH1_Lh_drP_total);
 EH2 = EH_local_170828_v3(eleLimit2, gasLimit2, EH2_Le, EH2_Lh, EH2_solarP, EH2_windP, CHP2_para, Boiler2_para, ES2_para, HS2_para, dev_L, dev_PV, dev_WT, EH2_solarP_rate, EH2_windP_rate, EH2_Le_drP_rate, EH2_Le_drP_total, EH2_Lh_drP_rate, EH2_Lh_drP_total);
 EH3 = EH_local_170828_v3(eleLimit3, gasLimit3, EH3_Le, EH3_Lh, EH3_solarP, EH3_windP, CHP3_para, Boiler3_para, ES3_para, HS3_para, dev_L, dev_PV, dev_WT, EH3_solarP_rate, EH3_windP_rate, EH3_Le_drP_rate, EH3_Le_drP_total, EH3_Lh_drP_rate, EH3_Lh_drP_total);

priceArray = elePrice;
num = 1;
for p1 = minMarketPrice : 0.01 : maxMarketPrice
    priceArray(1) = p1;
    [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice1, 1);
    demand1(num) = x(1);
    [x,~,~,~,~] = EH2.handlePrice(priceArray, gasPrice1, 1);
    demand2(num) = x(1);
    [x,~,~,~,~] = EH1.handlePrice(priceArray, gasPrice3, 1);
    demand3(num) = x(1);
    num = num + 1 ;
end
price = minMarketPrice : 0.01 : maxMarketPrice;
hold on
% plot(price,demand1,'LineWidth',1.5);
% plot(price,demand2,'LineWidth',1.5); 
% plot(price,demand3,'LineWidth',1.5);
% xlabel('电价')

plot(gasPrice1./price,demand1,'LineWidth',1.5);
plot(gasPrice1./price,demand2,'LineWidth',1.5); 
plot(gasPrice3./price,demand3,'LineWidth',1.5);
xlabel('气价/电价')
legend('IES1','IES2','IES3')
ylabel('需求')
