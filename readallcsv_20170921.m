clc
clear

%{
%读取负荷数据
file=dir('E:\JasonPC\Study\GridLAB-D\03-tool\PowerMatcher\多能源优化\博耳电力数据（282个）\*.csv');
period = 60/15; %优化周期是多少分钟，修改分母

count = 0; %count<n，只计数有效站点
% loadname={'initial','initial'};

for n=1:length(file)
    %     temp=dlmread(['E:\new\',file(n).name],' ',0,1);
    %     eval([file(n).name(1:end-4),'=temp;'])
    
    %打开csv文件
    fid = fopen(['E:\JasonPC\Study\GridLAB-D\03-tool\PowerMatcher\多能源优化\博耳电力数据（282个）\',file(n).name]);
    %读取表头 数据返回为cell类型 调用格式title{1}
    title = textscan(fid, '%s %s %s %s %s',1,'delimiter', ',');
    %读取数据 返回为cell类型
    data = textscan(fid, '%s %s %s %s %f','delimiter', ','); %d32是整数，s是字符串，f是浮点
    fclose(fid);
    
    power = data{1,5};
    if mod(length(power),(24*period))==0 && ~isempty(power) && min(power)>=0
        count = count + 1;
        loadName{count} = file(n).name;
        
        con = length(power) / (24*period);
        if con == 1
            loadValue(:,count) = power;
        else
            loadValue(1:(24*period) , count) = zeros(24*period , 1);
            for i = 1:length(power)
                loadValue(ceil(i/con) ,count) = loadValue(ceil(i/con) ,count) + power(i)/con;
            end
        end
        
        % 计算最高负荷、最低负荷所在的时间
        [peakPower, peakTime] = max( loadValue(:,count) );
%         maxPowerList(count,1) = maxPower;
        list_peakTime(count, 1) = (peakTime-1)/period;
%         if maxTimeList(count, 1) >= 8 && maxTimeList(count, 1) < 12 || maxTimeList(count, 1) >= 17 && maxTimeList(count, 1) < 21 || maxTimeList(count, 1) >= 0 && maxTimeList(count, 1) < 8
%             maxFlagList(count, 1) = 2;
%         else
%             maxFlagList(count, 1) = 1;
%         end
        [valleyPower, valleyTime] = min( loadValue(:,count) );
        list_valleyTime(count, 1) = (valleyTime-1)/period;
        
        
        % 负荷率：平均负荷与最高负荷的比率。
        avgPower = sum( loadValue(:,count) ) / length( loadValue(:,count) );
        list_avgRate(count, 1) = avgPower / peakPower;
                
        % 最小负荷率：报告期最低负荷与最高负荷的比率。
        list_valleyRate(count, 1) = valleyPower / peakPower;
        
        % 峰谷差：最高负荷与最低负荷之差。
        % 峰谷差率：峰谷差与最高负荷的比率。
        list_diffRate(count, 1) = (peakPower - valleyPower) / peakPower;
        
        %选取有峰谷差的电负荷
        %要求：最高负荷在峰时，最低负荷在谷时，三个比率适中
%         if ( list_peakTime(count, 1) >= 8 && list_peakTime(count, 1) < 12 || list_peakTime(count, 1) >= 17 && list_peakTime(count, 1) < 21 )...
%                 && ( list_valleyTime(count, 1) >= 0 && list_valleyTime(count, 1) < 8 )...
%                 && ( list_valleyRate(count, 1) >= 0.4 && list_valleyRate(count, 1) <= 0.6 )
%             list_flag(count, 1) = 11;
%         else
%             list_flag(count, 1) = 0;        
%         end
        
        %选取峰谷差很小的热负荷
        %要求：最高负荷在峰时，最低负荷在谷时，三个比率适中
        if ( list_peakTime(count, 1) >= 8 && list_peakTime(count, 1) < 12 || list_peakTime(count, 1) >= 17 && list_peakTime(count, 1) < 21 )...
                && ( list_valleyTime(count, 1) >= 0 && list_valleyTime(count, 1) < 8 )...
                && ( list_valleyRate(count, 1) >= 0.8 )
            list_flag(count, 1) = 11;
        else
            list_flag(count, 1) = 0;        
        end
        
    else
        disp([file(n).name,' 舍弃该数据'])
    end
    
    
end

% 之前存的1h间隔的数据
% save tmp1.mat loadName
% save tmp2.mat loadValue
%20171229 改为存储30min数据
% save data_loadName_30min.mat loadName
% save data_loadValue_30min.mat loadValue
%20180201 改为存储15min数据
save data_loadName_15min.mat loadName
save data_loadValue_15min.mat loadValue
%}


%读取可再生能源的数据
file=dir('E:\JasonPC\Study\GridLAB-D\03-tool\PowerMatcher\Data\刘东_天气-功率转换功率\刘东团队_天气-功率转换功率\Penn_State_PA_result\*.txt');
period = 60 / 15; %优化周期是多少分钟，修改分母

count = 0; %count<n，只计数有效站点
% loadname={'initial','initial'};

for n=1:length(file)
    %     temp=dlmread(['E:\new\',file(n).name],' ',0,1);
    %     eval([file(n).name(1:end-4),'=temp;'])
    
    %打开csv文件
    fid = fopen(['E:\JasonPC\Study\GridLAB-D\03-tool\PowerMatcher\Data\刘东_天气-功率转换功率\刘东团队_天气-功率转换功率\Penn_State_PA_result\',file(n).name]);
    %读取表头 数据返回为cell类型 调用格式title{1}
    title = textscan(fid, '%s %s %s %s %s %s %s %s',1,'delimiter', '\t');
    %读取数据 返回为cell类型
    data = textscan(fid, '%s %s %f %f %f %f %f %f','delimiter', '\t'); %d32是整数，s是字符串，f是浮点
    fclose(fid);
    
    solar = data{1,6};
    wind = data{1,5};
    
    
    
    if mod(length(solar),(24*period))==0 && ~isempty(solar) && min(solar)>=0 && max(solar)<2e5 ... %谨慎使用max条件
        && mod(length(wind),(24*period))==0 && ~isempty(wind) && min(wind)>=0
    
        count = count + 1;
        renewableName{count} = file(n).name;
        
        con = length(solar) / (24*period);
        if con == 1
            solarValue(:,count) = solar;
        else
            solarValue(1:(24*period) , count) = zeros(24*period , 1);
            for i = 1:length(solar)
                solarValue(ceil(i/con) ,count) = solarValue(ceil(i/con) ,count) + solar(i)/con;
            end
        end
        
     
        con = length(wind) / (24*period);
        if con == 1
            windValue(:,count) = wind;
        else
            windValue(1:(24*period) , count) = zeros(24*period , 1);
            for i = 1:length(wind)
                windValue(ceil(i/con) ,count) = windValue(ceil(i/con) ,count) + wind(i)/con;
            end
        end

    else
        disp([file(n).name,' 舍弃该风电与光伏数据'])
    end

end

% 之前存的1h间隔的数据
% save renewableName.mat renewableName
% save solarValue.mat solarValue
% save windValue.mat windValue
%20171229 改为存储30min数据
% save renewableName_30min.mat renewableName
% save solarValue_30min.mat solarValue
% save windValue_30min.mat windValue
%20180201 改为存储15min数据
save renewableName_15min.mat renewableName
save solarValue_15min.mat solarValue
save windValue_15min.mat windValue