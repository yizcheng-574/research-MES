figure;
IES_no = 2;
P_charge = zeros(24, 3);%1-3列分别为P2最优解，P2transform; P1最优解，
P_discharge = zeros(24, 3);
P_wind = zeros(24, 3);
green = [59 129 113]/255;
for isMILP = 0 : 1
    ies_no = IES_no + isMILP * IESNUMBER;
    %P2和P1的最优风功率
    P_wind(:, 2 * isMILP + 1) = - result_Ele_loss(:, ies_no) -  result_CHP_power(:, ies_no)...
        - result_ES_discharge(:, ies_no) + result_ES_charge(:, ies_no) + result_EH_Le(:, ies_no) + result_EH_Edr(:, ies_no);
    %P2和P1最优充电功率
    P_charge(: , 2 * isMILP + 1) = result_ES_charge(:, ies_no);
    P_discharge(:, 2 * isMILP + 1) = result_ES_discharge(:, ies_no);
end
%MPPT风功率曲线
P_renewable = result_EH_windP(:, IES_no) + result_EH_solarP(:, IES_no);

delta_S = P_charge * 0.9 - P_discharge / 0.9;
%P2修正后的充放电
for i = 1 : 24
    if delta_S(i, 1) > 0
        P_charge(i, 2) = delta_S(i, 1) / 0.9;
    elseif delta_S(i, 1) < 0
        P_discharge(i, 2) = - delta_S(i, 1) * 0.9;
    end
end
%P2修正后的弃风
P_wind(:, 2) = - result_Ele_loss(:, IES_no) -  result_CHP_power(:, IES_no)...
   - P_discharge(:, 2) + P_charge(:, 2) + result_EH_Le(:, IES_no) + result_EH_Edr(:, IES_no);
figure; hold on;
w = 2;

yyaxis left;
H1 = bar(t, -P_charge);
H2 = bar(t, P_discharge);
ylabel('EES power(kW)');
H1(1).EdgeColor = brown;
H1(1).FaceColor = H1(1).EdgeColor;
H1(2).EdgeColor = gray;
H1(2).FaceColor = H1(2).EdgeColor;
H1(3).EdgeColor = green;
H1(3).FaceColor = H1(3).EdgeColor;
H2(1).EdgeColor = brown;
H2(1).FaceColor = H1(1).EdgeColor;
H2(2).EdgeColor = gray;
H2(2).FaceColor = H1(2).EdgeColor;
H2(3).EdgeColor = green;
H2(3).FaceColor = H2(3).EdgeColor;

yyaxis right;
t1 = 1:25;
H3 = stairs(t1, appendStairArray(P_renewable) / 1000, 'LineWidth', w);
H4 = stairs(t1, appendStairArray(P_wind) / 1000, 'LineWidth', w);
set(H3, 'Color', gold, 'LineWidth', 3);
set(H4(1), 'Color', brown, 'LineStyle','--', 'Marker', '.', 'MarkerSize', 13);
set(H4(2), 'Color', gray, 'LineStyle','--', 'Marker', '.', 'MarkerSize', 13);
set(H4(3), 'Color', green, 'LineStyle','--','Marker', 'd', 'MarkerSize', 13);

ylabel('renewable power(kW)');

le = legend([H3, H4(1), H4(2), H4(3), H1(1), H1(2), H1(3)], ...
    'MPPT output', 'optimal power of (P2)','after transformation', 'optimal power of (P1)',...
    'EES power of (P2)', 'EES power after transformation', 'EES power of (P1)', ...
     'Location','northoutside','Orientation','vertical');
set(le,'Box','off');
legend('show');
xlabel(sprintf('MES%d', IES_no));
xlim([0, 24 * period + 1]);
xticks(0:(24 * period / 4) : 24 * period);
xticklabels({ '0:00','6:00','12:00','18:00','24:00' });
set(gcf,'Position',[0 0 500 500]);
