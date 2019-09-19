function [accomodation_rate ] = calculate_accomodation_rate(result_Ele, gridClearDemand, EH_res_total, res_total, result_EH_windP, result_EH_solarP )
    unaccomodated_res = EH_res_total;
    unaccomodated_res_total = 0;

    for tmp_t = 1 : 24
        for tmp_ies = 1 : 3
            if result_Ele(tmp_t, tmp_ies) < 0
                unaccomodated_res(tmp_t) = unaccomodated_res(tmp_t) +...
                    min(-result_Ele(tmp_t, tmp_ies), result_EH_windP(tmp_t, tmp_ies) + result_EH_solarP(tmp_t, tmp_ies));
            end
        end
        if gridClearDemand(tmp_t) < 0
            unaccomodated_res_total = unaccomodated_res_total + ...
                    min(-gridClearDemand(tmp_t), unaccomodated_res(tmp_t) + EH_res_total(tmp_t));
        end
    end
    accomodation_rate = unaccomodated_res_total / res_total;
end

