function [ result ] = appendStairArray(array)%为stairs作图的数组增加最后一列
    [row, col] = size(array);
    if row == 1 %行向量
        result = [array, array(end)];
    elseif col == 1 %列向量
        result = [array; array(end)];
    else %矩阵
        result = [array; array(end, :)];
    end
end