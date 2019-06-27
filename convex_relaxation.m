% IES_no = 2;
% P_charge = zeros(24, 3);%1-3列分别为P2最优解，P2transform; P1最优解，
% P_discharge = zeros(24, 3);
% P_wind = zeros(24, 3);
% green = [59 129 113]/255;
% for isMILP = 0 : 1
%     ies_no = IES_no + isMILP * IESNUMBER;
%     %P2和P1的最优风功率
%     P_wind(:, 2 * isMILP + 1) = - result_Ele_loss(:, ies_no) -  result_CHP_power(:, ies_no)...
%         - result_ES_discharge(:, ies_no) + result_ES_charge(:, ies_no) + result_EH_Le(:, ies_no) + result_EH_Edr(:, ies_no);
%     %P2和P1最优充电功率
%     P_charge(: , 2 * isMILP + 1) = result_ES_charge(:, ies_no);
%     P_discharge(:, 2 * isMILP + 1) = result_ES_discharge(:, ies_no);
% end
% %MPPT风功率曲线
% P_renewable = result_EH_windP(:, IES_no) + result_EH_solarP(:, IES_no);
% 
% delta_S = P_charge * 0.9 - P_discharge / 0.9;
% %P2修正后的充放电
% for i = 1 : 24
%     if delta_S(i, 1) > 0
%         P_charge(i, 2) = delta_S(i, 1) / 0.9;
%     elseif delta_S(i, 1) < 0
%         P_discharge(i, 2) = - delta_S(i, 1) * 0.9;
%     end
% end
% %P2修正后的弃风
% P_wind(:, 2) = - result_Ele_loss(:, IES_no) -  result_CHP_power(:, IES_no)...
%    - P_discharge(:, 2) + P_charge(:, 2) + result_EH_Le(:, IES_no) + result_EH_Edr(:, IES_no);
colorP1 = ColorHex(['F5','48','45'])/255;
% colorP2 = ColorHex(['37','76','F5'])/255;
colorP2 = dodgerblue;
colorP2EEC = ColorHex(['51','CA','52'])/255;

figure; hold on;
w = 2;
yyaxis left;
H1 = bar(t, -[P_charge(:,3), P_charge(:,1), P_charge(:,2)]);
H2 = bar(t, [P_discharge(:,3), P_discharge(:,1), P_discharge(:,2)]);
ylabel('EES power(kW)');
H1(1).EdgeColor = colorP1;
H1(1).FaceColor = H1(1).EdgeColor;
H1(2).EdgeColor = colorP2;
H1(2).FaceColor = H1(2).EdgeColor;
H1(3).EdgeColor = colorP2EEC;
H1(3).FaceColor = H1(3).EdgeColor;
H2(1).EdgeColor = colorP1;
H2(1).FaceColor = H1(1).EdgeColor;
H2(2).EdgeColor = colorP2;
H2(2).FaceColor = H1(2).EdgeColor;
H2(3).EdgeColor = colorP2EEC;
H2(3).FaceColor = H2(3).EdgeColor;
ylim([-200, 200]);
yticks(-200:100:200);
yticklabels({'-200', '-100', '0', '100', '200'});
yyaxis right;
t1 = 1:25;
H4 = stairs(t1, appendStairArray(P_renewable - P_wind(:, 1)), 'LineWidth', w);
set(H4, 'Color', 'k', 'LineStyle','-', 'Marker', '.', 'MarkerSize', 13);
ylim([-800, 800]);
yticks(-800:800:800);
yticklabels({'-800', '0', '800'});
ylabel('power curtailment(kW)');

le = legend([H4, H1(1), H1(2), H1(3)], ...
    'curtailed power of (P2)',...
    '(P1)', '(P2)', '(P2+EEC)', ...
     'Location','northoutside','Orientation','vertical');
set(le,'Box','off');
legend('show');
xlabel(sprintf('MES%d', IES_no));
xlim([0, 24 * period + 1]);
xticks(0:(24 * period / 4) : 24 * period);
xticklabels({ '0:00','6:00','12:00','18:00','24:00' });
set(gcf,'Position',[0 0 500 500]);
