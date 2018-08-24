
function clearPrice = clearing(demandArray, targetDemand) %一般targetDemand=0

global minMarketPrice maxMarketPrice step

leftIx = 1;
rightIx = length(demandArray);
%demandArray是所有投标的聚合，且总是递减的

%First test for a few special cases
if targetDemand > demandArray(leftIx)
    %If the target is higher than the maximum of the bid, return the minimum price
    clearPrice = minMarketPrice;
    return % uses return to handle the special case
elseif targetDemand < demandArray(rightIx)
    % If the target is lower than the minimum of the bid, return the maximum price
    clearPrice = maxMarketPrice;
    return
elseif targetDemand == demandArray(leftIx)  %如果等于左端点，则右端点移到左端点
    rightIx = leftIx; % 跳过二分法那个while环节
elseif targetDemand == demandArray(rightIx)
    leftIx = rightIx; % 跳过二分法那个while环节
else % demand is between the limits of this bid, which can not be flat at this point
    % Go on while there is at least 1 point between the left and right index
    % 二分法
    while rightIx - leftIx > 1
        % Determine the middle between the 2 boundaries
        middleIx = fix((leftIx + rightIx) / 2);
        middleDemand = demandArray(middleIx);
        
        if targetDemand==middleDemand
            % A point with the target demand is found, select this point
            leftIx = middleIx;
            rightIx = middleIx;
        elseif middleDemand > targetDemand
            % If the middle demand is bigger than the target demand, we set the left to the middle
            leftIx = middleIx;
        else % middleDemand < targetDemand
            % If the middle demand is smaller than the target demand, we set the right to the middle
            rightIx = middleIx;
        end
    end
end

% If the left or right point matches the targetDemand, expand the range
% 扩展等于target的范围
while leftIx > 1 && targetDemand == demandArray(leftIx - 1)
    leftIx = leftIx - 1;
end
while rightIx < length(demandArray) && targetDemand == demandArray(rightIx + 1)
    rightIx = rightIx + 1;
end
%此时的左端点到右端点之间所有的点都满足target，或左右相邻都不等于target但包含target

% return interpolate(leftIx, rightIx, targetDemand);
% 将interpolate方法直接加入到本方法中

% 三目运算符 x? y:z
% x是一个boolean类型,若x为true,结果显示y,若x为false,则结果显示z.
if rightIx == 1
    leftPrice = minMarketPrice;
else
    leftPrice = minMarketPrice + (leftIx-1) * step;
end

if leftIx == length(demandArray)
    rightPrice = maxMarketPrice;
else
    rightPrice = minMarketPrice + (rightIx-1) * step;
end

leftDemand = demandArray(leftIx);
rightDemand = demandArray(rightIx);

if leftDemand == rightDemand
    demandFactor = 0.5;
else
    demandFactor = (leftDemand - targetDemand) / (leftDemand - rightDemand);
end

clearPrice = leftPrice + (rightPrice - leftPrice) * demandFactor;

end

% test
% targetDemand=0;
% demandArray=[100; 90; 80; 70; 60; 50; 40; 30; 20; 10; 9];
% demandArray=[100; 90; 80; 70; 60; 50; 40; 30; 20; 20; -10];
% demandArray=[100; 90; 80; 70; 60; 50; 40; 30; 0; 0; -10];
% demandArray=[100; 90; 80; 70; 60; 50; 40; 30; 0; 0; 0];
% demandArray=[-100; -290; -380; -470; -560; -650; -740; -830; -830; -830; -830];
% demandArray=[0; 0; 0; -470; -560; -650; -740; -830; -830; -830; -830];
% clearPrice = clearing(demandArray, targetDemand)
