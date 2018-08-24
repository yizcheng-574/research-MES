% 预测业务
% 20180219 删除了seedNumber参数
% 风光的标准差分开

%{
function [Le_result, Lh_result] = predict(number, time) % EH的编号，第几小时的预测，time=0是日前，1-24是日内
    switch number
        case 1
            Le_base = [115
                105
                102
                103
                105
                110
                120
                125
                122
                116
                115
                110
                116
                119
                120
                128
                129
                130
                132
                136
                130
                120
                110
                108]; %日前预测负荷
            
            Lh_base = [110
                109
                108
                107
                106
                105
                106
                107
                108
                107
                106
                105
                103
                102
                102
                104
                106
                110
                111
                113
                115
                114
                112
                110].*1;
        otherwise
            Le_base = zeros(24,1);
            Lh_base = zeros(24,1);
    end
    
    if time == 0
        Le_result = Le_base;
        Lh_result = Lh_base;
    else
        Le_error = randn([(24+1-time),1])*1; %更新
        Lh_error = randn([(24+1-time),1])*1; %更新
        if time ~= 1
            Le_result(1,:)=[]; %更新
            Lh_result(1,:)=[]; %更新
        end
        Le_result = Le_result + Le_error;
        Lh_result = Lh_result + Lh_error;
    end

end
%}

function [Le_result, Lh_result, solarP_result, windP_result] = predict(Le, Lh, solarP, windP, t_current, dev_L, dev_PV, dev_WT, solarP_rate, windP_rate) % t=1-24都是日内，没有日前
%     global period

    %rand 生成均匀分布的伪随机数 分布在（0~1）之间
    %randn 生成标准正态分布的伪随机数 （均值为0，方差为1）
%     randn('seed', t_current);
    
    %电、热负荷
%     Le_error = zeros(24*period,1);
%     Lh_error = zeros(24*period,1);
%     Le_error(t_current:24*period) = randn([(24*period-t_current+1),1]) .* Le(t_current:24*period) * dev_L; %预测误差
%     Lh_error(t_current:24*period) = randn([(24*period-t_current+1),1]) .* Lh(t_current:24*period) * dev_L;
%     Le_result = Le + Le_error; %更新
%     Lh_result = Lh + Lh_error;
    %改为只计当前时刻的预测误差，否则累计的误差过大；频繁重复预测也不科学
    Le_error = randn() * Le(t_current) * dev_L; %预测误差
    Lh_error = randn() * Lh(t_current) * dev_L;
    Le_result = Le;
    Lh_result = Lh;
    Le_result(t_current) = Le_result(t_current) + Le_error; %更新
    Lh_result(t_current) = Lh_result(t_current) + Lh_error;
    
    %风、光
    %有几点和负荷不一样
    %一是：如果光是零，那么一定是零，没有预测误差
    %二是：风、光接近零的时候，不要随机变为负数
%     solarP_error = zeros(24*period,1);
%     windP_error = zeros(24*period,1);
%     solarP_error(t_current:24*period) = randn([(24*period-t_current+1),1]) .* solarP(t_current:24*period) * dev_RES; %预测误差
%     windP_error(t_current:24*period) = randn([(24*period-t_current+1),1]) .* windP(t_current:24*period) * dev_RES;
%     %修正
%     for i = t_current : 24*period
%         if solarP(i) == 0
%             solarP_error(i) = 0;
%         end
%         if solarP(i) + solarP_error(i) < 0
%             solarP_error(i) = - solarP(i);
%         end
%         if windP(i) == 0
%             windP_error(i) = 0;
%         end
%         if windP(i) + windP_error(i) < 0
%             windP_error(i) = - windP(i);
%         end
%     end   
%     solarP_result = solarP + solarP_error; %更新
%     windP_result = windP + windP_error;
    %改为只计当前时刻的预测误差，否则累计的误差过大；频繁重复预测也不科学
    solarP_error = randn() * solarP(t_current) * dev_PV; %预测误差，solarP(t_current)=0的时候自然为零
    windP_error = randn() * windP(t_current) * dev_WT; %预测误差
    
    solarP_result = solarP;
    windP_result = windP;
    
    solarP_result(t_current) = solarP_result(t_current) + solarP_error; %更新
    if solarP_result(t_current) < 0
        solarP_result(t_current) = 0;
    elseif solarP_result(t_current) > solarP_rate
        solarP_result(t_current) = solarP_rate;
    end
    
    windP_result(t_current) = windP_result(t_current) + windP_error;
    if windP_result(t_current) < 0
        windP_result(t_current) = 0;
    elseif windP_result(t_current) > windP_rate
        windP_result(t_current) = windP_rate;
    end

end
