clc;
clear;
close all

StartYear = 1961; % 数据起始年份
EndYear = 2022; % 数据结束年份

Year = StartYear:EndYear;
ID = ['E:\北大兼职\臭氧高温复合\extreme_pre\', 'CHM_PRE_0.1dg_19612022.nc'];
lat = ncread(ID, 'latitude');
lon = ncread(ID, 'longitude');
year = Year - StartYear + 1;

% 提取所需地区的坐标
lat = lat(16:86);
lon = lon(366:466);
stlo = 366; stla = 16; stpre = 1; 
start = [stlo, stla, stpre]; % 创建start变量，每一维开始的变量
locount = 101; lacount = 71; ticount = 22645;
count = [locount, lacount, ticount]; % 创建count变量，从每一维的start开始读取的总数目
stride = [1, 1, 1]; % 设置读取的步长
pre_spei_stride1 = ncread(ID, 'pre', start, count, stride);
pre = double(pre_spei_stride1);

cdL_NMG=shaperead("E:\北大兼职\臭氧高温复合\extreme_pre\广东省\广东省.shp");
cx_NMG=[cdL_NMG(:).X]; cy_NMG=[cdL_NMG(:).Y];

for i = 1:101
    for j = 1:71
        in = inpolygon(lon(i),lat(j),cx_NMG,cy_NMG);
        if in == 0 % 在边界范围外
            pre(i,j,:) = NaN;
        end
    end
end

% 获取每年的天数
days_in_year = @(y) 365 + (mod(y, 4) == 0 & (mod(y, 100) ~= 0 | mod(y, 400) == 0));

% 初始化极端降水指数
PRCPTOT_seasonal = zeros(locount, lacount, 4, length(Year));
R10_seasonal = zeros(locount, lacount, 4, length(Year));
R20_seasonal = zeros(locount, lacount, 4, length(Year));

% 定义季节月份
seasons = {3:5, 6:8, 9:11, [12, 1, 2]};

% 遍历每年
for y = 1:length(Year)
    % 计算当前年之前的总天数
    days_before = sum(arrayfun(days_in_year, Year(1:y-1)));
    
    % 获取当前年的天数
    num_days = days_in_year(Year(y));
    
    % 获取当前年的降水数据
    pre_year = pre(:, :, days_before + 1 : days_before + num_days);
    
    % 按季节遍历
    for s = 1:4
        % 获取当前季节的月份
        months = seasons{s};
        
        % 初始化季节降水数据
        pre_season = [];
        
        % 为冬季处理不同年份的数据
        if s == 4
            for m = months
                if m == 12
                    pre_season = cat(3, pre_season, pre_year(:, :, sum(eomday(Year(y), 1:11)) + 1:sum(eomday(Year(y), 1:12))));
                elseif m == 1 || m == 2
                    % 处理下一年的1月和2月，检查是否有下一年
                    if y < length(Year)
                        days_next_year = days_in_year(Year(y+1));
                        days_before_next_year = sum(arrayfun(days_in_year, Year(1:y)));
                        pre_next_year = pre(:, :, days_before_next_year + 1 : days_before_next_year + days_next_year);
                        start_day = sum(eomday(Year(y+1), 1:m-1)) + 1;
                        end_day = start_day + eomday(Year(y+1), m) - 1;
                        pre_season = cat(3, pre_season, pre_next_year(:, :, start_day:end_day));
                    end
                end
            end
        else
            for m = months
                start_day = sum(eomday(Year(y), 1:m-1)) + 1;
                end_day = start_day + eomday(Year(y), m) - 1;
                pre_season = cat(3, pre_season, pre_year(:, :, start_day:end_day));
            end
        end
        
        % 计算当前季节的极端降水指数
        for i = 1:locount
            for j = 1:lacount
                PRCPTOT_seasonal(i, j, s, y) = sum(pre_season(i, j, pre_season(i, j, :) > 1), 'all');
                R10_seasonal(i, j, s, y) = sum(pre_season(i, j, :) >= 10, 'all');
                R20_seasonal(i, j, s, y) = sum(pre_season(i, j, :) >= 20, 'all');
            end
        end
    end
end

for i = 1:101
    for j = 1:71
        for s = 1:4
            in = inpolygon(lon(i), lat(j), cx_NMG, cy_NMG);
            if in == 0 % 在边界范围外
                PRCPTOT_seasonal(i, j, s, :) = NaN;
                R10_seasonal(i, j, s, :) = NaN;
                R20_seasonal(i, j, s, :) = NaN;
            end
        end
    end
end

% 创建NetCDF文件并写入数据
ncid = netcdf.create('E:\北大兼职\臭氧高温复合\extreme_pre\extreme_precipitation_indices_seasonal.nc', 'NETCDF4');

% 定义维度
dimid_lon = netcdf.defDim(ncid, 'lon', locount);
dimid_lat = netcdf.defDim(ncid, 'lat', lacount);
dimid_season = netcdf.defDim(ncid, 'season', 4);
dimid_year = netcdf.defDim(ncid, 'year', length(Year));

% 定义变量
varid_lon = netcdf.defVar(ncid, 'lon', 'double', dimid_lon);
varid_lat = netcdf.defVar(ncid, 'lat', 'double', dimid_lat);
varid_year = netcdf.defVar(ncid, 'year', 'double', dimid_year);
varid_prcptot = netcdf.defVar(ncid, 'PRCPTOT', 'double', [dimid_lon, dimid_lat, dimid_season, dimid_year]);
varid_r10 = netcdf.defVar(ncid, 'R10', 'double', [dimid_lon, dimid_lat, dimid_season, dimid_year]);
varid_r20 = netcdf.defVar(ncid, 'R20', 'double', [dimid_lon, dimid_lat, dimid_season, dimid_year]);

% 结束定义模式
netcdf.endDef(ncid);

% 写入数据
netcdf.putVar(ncid, varid_lon, lon);
netcdf.putVar(ncid, varid_lat, lat);
netcdf.putVar(ncid, varid_year, Year);
netcdf.putVar(ncid, varid_prcptot, PRCPTOT_seasonal);
netcdf.putVar(ncid, varid_r10, R10_seasonal);
netcdf.putVar(ncid, varid_r20, R20_seasonal);

% 关闭NetCDF文件
netcdf.close(ncid);
