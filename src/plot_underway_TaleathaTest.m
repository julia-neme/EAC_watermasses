%Plot ship underwater data 

path(path,genpath('./Utilities'));
path(path,genpath('./Utilities/export_fig'));

dir='../../underway'
file='in2022_v06uwy.nc';
% % % % % % dir='../../../../IN2015_V02/underway'
% % % % % % file='in2015_v02uwy.nc';

%% LOAD DATA Coast and bathy
load BATHY/eaccoast.dat
file_bathy='BATHY/bathy_dbdb2_v30_AUSTRALIA.nc';
lon_bat=ncread(file_bathy,'lon');
lat_bat=ncread(file_bathy,'lat');
height_bat=ncread(file_bathy,'height');

%% LOAD DATA ADCP file
path_file=[dir, '/', file];
ncdisp(path_file) % ncdump(file)

ncid = netcdf.open (path_file);
[numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid);

% start time = 13211755 seconds since 2015/01/01
% Global Attributes:
%            Conventions          = 'COARDS'
%            Ship                 = 'Investigator'
%            Voyage               = 'in2022_v06'
%            Epoch                = '16762310 seconds since 2022-01-01 00:00:00 00'
%            SampleInterval       = 5
%            MaxInterpolationSpan = 24
%            SamplingUnits        = 'Sampling interval is in seconds, max. interpolation span is in samples'

start_days=16762310/(60*60*24);
end_days = start_days + (5 * (len_time-1))/(60*60*24);
start_time = datenum(2022,1,1,0,0,0)+start_days; % change to 2022?
end_time = datenum(2022,1,1,0,0,0)+end_days; % change to 2022?
len_time_temp = size(lat); len_time = len_time_temp(1);
time = start_time:datenum(0,0,0,0,0,5):end_time;

% lon lat depth
lon_vinfo = ncinfo(path_file,'longitude')
lat_vinfo = ncinfo(path_file,'latitude')
depth_vinfo = ncinfo(path_file,'depth')
% time_vinfo = ncinfo(path_file,'time')
lon = ncread(path_file,'longitude');   lon(lon>1.e+37 | lon==0)=NaN;
lat = ncread(path_file,'latitude');   lat(lat>1.e+37 | lat==0)=NaN;
depth = double(ncread(path_file,'depth')); depth(depth>1.e+37)=NaN;
% time = ncread(path_file,'time')+datenum(2015,01,01);
% time = ncread(path_file,'time');



% waterFlow
waterFlow_vinfo = ncinfo(path_file,'waterFlow')
waterFlow_vinfo.Attributes(1,1).Value   % long name
waterFlow_vinfo.Attributes(1,2).Value   % units
waterFlow = ncread(path_file,'waterFlow'); 

% temp
temp_vinfo = ncinfo(path_file,'waterTemp')
temp_vinfo.Attributes(1,1).Value   % long name
temp_vinfo.Attributes(1,2).Value   % units
temp = ncread(path_file,'waterTemp');   temp(temp>1.e+37)=NaN;    temp(waterFlow<2 | isnan(waterFlow))=NaN;  

% salinity
sal_vinfo = ncinfo(path_file,'salinity')
sal_vinfo.Attributes(1,1).Value   % long name
sal_vinfo.Attributes(1,2).Value   % units
sal = ncread(path_file,'salinity');   sal(sal>1.e+37)=NaN;    sal(waterFlow<2 | isnan(waterFlow))=NaN;  
sal(sal<35.4)=NaN;  

% fluorescence
flr_vinfo = ncinfo(path_file,'fluorescence')
flr_vinfo.Attributes(1,1).Value   % long name
flr_vinfo.Attributes(1,2).Value   % units
flr = ncread(path_file,'fluorescence');   flr(flr>1.e+37)=NaN;    flr(waterFlow<2 | isnan(waterFlow))=NaN;  

% portTrueWindDir
WindDir_port_vinfo = ncinfo(path_file,'portTrueWindDir')
WindDir_port_vinfo.Attributes(1,1).Value   % long name
WindDir_port_vinfo.Attributes(1,2).Value   % units
WindDir_port = ncread(path_file,'portTrueWindDir');   WindDir_port(WindDir_port>1.e+37)=NaN; 
% rawPortTrueWindSpeed
WindSpeed_port_vinfo = ncinfo(path_file,'rawPortTrueWindSpeed')
WindSpeed_port_vinfo.Attributes(1,1).Value   % long name
WindSpeed_port_vinfo.Attributes(1,2).Value   % units
WindSpeed_port = ncread(path_file,'rawPortTrueWindSpeed');   WindSpeed_port(WindSpeed_port>1.e+37)=NaN; 

% stbdTrueWindDir
WindDir_stbd_vinfo = ncinfo(path_file,'stbdTrueWindDir')
WindDir_stbd_vinfo.Attributes(1,1).Value   % long name
WindDir_stbd_vinfo.Attributes(1,2).Value   % units
WindDir_stbd = ncread(path_file,'stbdTrueWindDir');   WindDir_stbd(WindDir_stbd>1.e+37)=NaN; 
% stbdTrueWindSpeed
WindSpeed_stbd_vinfo = ncinfo(path_file,'stbdTrueWindSpeed')
WindSpeed_stbd_vinfo.Attributes(1,1).Value   % long name
WindSpeed_stbd_vinfo.Attributes(1,2).Value   % units
WindSpeed_stbd = ncread(path_file,'stbdTrueWindSpeed');   WindSpeed_stbd(WindSpeed_stbd>1.e+37)=NaN;   

%% Plot temperature , & salinity time series
figure
subplot(3,1,1)
plot(time,temp,'.b')
datetick('x','hh/dd/yyyy','keeplimits')
ylabel('temperature (^oC)','fontsize',12)
%xlabel('time (day/hour)','Fontsize',12)
set(gca,'xticklabel',[])
set(gca,'Fontsize',12)

subplot(3,1,2)
plot(time,sal,'.r')
datetick('x','hh/dd/yyyy','keeplimits')
ylabel('Salinity ()','fontsize',12)
%xlabel('time (day/hour)','Fontsize',12)
%set(gca,'xticklabel',[])
set(gca,'Fontsize',12)

subplot(3,1,3)
plot(time,flr,'.g')
datetick('x','hh/dd/yyyy','keeplimits')
ylabel('Fluorescence ()','fontsize',12)
%xlabel('time (day/hour)','Fontsize',12)
%set(gca,'xticklabel',[])
set(gca,'Fontsize',12)


%% Plot wind time series
figure
set(gcf, 'position',[4   4   1200   700],'PaperPositionMode','auto' , 'Color', 'w')
s3=subplot(2,1,1)
plot(time,WindSpeed_port)
hold on
plot(time,WindSpeed_stbd,'r')
set(gca,'Fontsize',12)
datetick('x','dd mmm')
xlabel('Date (2022)')
title('Wind speed (underway)')
ylabel('KNOTS')
ylim([0 50])
legend('port','starboard')
grid on

s2=subplot(2,1,2)
plot(time,WindDir_port)
hold on
plot(time,WindDir_stbd,'r')
set(gca,'Fontsize',12)
datetick('x','dd mmm')
xlabel('Date (2022)')
title('Wind direction (underway)')
ylabel('KNOTS')
ylim([0 380])
legend('port','starboard')
grid on

export_fig(['../FIGURES/underway_Wind_'  datestr(time(length(temp)),'dd_HHMM') '.png'],'-opengl','-m2','-q101')





%% plot temp sal fluo
scale=2;
pos=[ 182         202        1319         504];
figure
set(gcf, 'position',pos,'PaperPositionMode','auto' , 'Color', 'w')
subplot(1,3,1)
scatter(lon(1:floor(length(lon)/500):end),lat(1:floor(length(lon)/500):end),40,temp(1:floor(length(lon)/500):end),'filled'); 
colorbar; caxis([20 24]) %caxis([21.5 24.5])
hold on; plot(eaccoast(:,1),eaccoast(:,2),'-k','LineWidth',2);
contour(lon_bat,lat_bat,height_bat',[-2000 -200],'color',[0.4 0.4 0.4]);
xlim([151 156]);
ylim([-33 -26.7])
Y=get(gca,'ylim');
set(gca,'dataaspectratio',[1 cos((Y(2)+Y(1))/2/180*pi) 1])
title(['Underway temperature ' file ],'Fontsize',14,'Interpreter','none')
xlabel('Longitude','Fontsize',14)
ylabel('Latitude','Fontsize',14)
set(gca,'Fontsize',12)
colormap(jet)

subplot(1,3,2)
scatter(lon(1:floor(length(lon)/500):end),lat(1:floor(length(lon)/500):end),40,sal(1:floor(length(lon)/500):end),'filled'); colorbar; caxis([35.55 35.75])
hold on; plot(eaccoast(:,1),eaccoast(:,2),'-k','LineWidth',2);
contour(lon_bat,lat_bat,height_bat',[-2000 -200],'color',[0.4 0.4 0.4]);
xlim([151 156]);
ylim([-33 -26.7])
Y=get(gca,'ylim');
set(gca,'dataaspectratio',[1 cos((Y(2)+Y(1))/2/180*pi) 1])
title(['Underway salinity ' file ],'Fontsize',14,'Interpreter','none')
xlabel('Longitude','Fontsize',14)
ylabel('Latitude','Fontsize',14)
set(gca,'Fontsize',12)
colormap(jet)

subplot(1,3,3)
scatter(lon(1:floor(length(lon)/500):end),lat(1:floor(length(lon)/500):end),40,flr(1:floor(length(lon)/500):end),'filled'); colorbar; caxis([1 3])
hold on; plot(eaccoast(:,1),eaccoast(:,2),'-k','LineWidth',2);
contour(lon_bat,lat_bat,height_bat',[-2000 -200],'color',[0.4 0.4 0.4]);
xlim([151 156]);
ylim([-33 -26.7])
Y=get(gca,'ylim');
set(gca,'dataaspectratio',[1 cos((Y(2)+Y(1))/2/180*pi) 1])
title(['Underway fluorescence ' file ],'Fontsize',14,'Interpreter','none')
xlabel('Longitude','Fontsize',14)
ylabel('Latitude','Fontsize',14)
set(gca,'Fontsize',12)
colormap(jet)

%export_fig(strcat('../plots/underway_06_23h.png'),'-painters','-m2','-q101')
%export_fig(strcat('../plots/underway_08.png'),'-opengl','-m2','-q101')






