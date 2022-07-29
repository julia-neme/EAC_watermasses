clear all
close all

science_path='X:\current\science\'
%load([science_path,'\DATA\Satellite\202105100833.mat'])
%filename_sst = [science_path,'\DATA\Satellite\SST\Brisbane2_4hr\2022072414.mat']
filename_sst = [science_path,'\DATA\Satellite\SST\Brisbane2_4hr\2022072422.mat']
%%filename_sst = [science_path,'\DATA\Satellite\SST\Brisbane2_4hr\2022071822.mat']
sst = load(filename_sst)
date_sst = filename_sst(54:63)

% addpath
addpath([science_path,'\Scripts\Utilities\cmocean\']);
addpath([science_path,'\Scripts\Utilities\altmany-export_fig-7720793']);

%%%%%%%% bathy
lon_bat=ncread('./BATHY/bathy_dbdb2_v30_AUSTRALIA.nc','lon') ;
lat_bat=ncread('./BATHY/bathy_dbdb2_v30_AUSTRALIA.nc','lat') ;
height_bat=ncread('./BATHY/bathy_dbdb2_v30_AUSTRALIA.nc','height') ;
load ./BATHY/eaccoast.dat
%%%%%%%%

%%%
%filename = ['https://thredds.aodn.org.au/thredds/dodsC/IMOS/OceanCurrent/GSLA/NRT00/2022/IMOS_OceanCurrent_HV_20220716T060000Z_GSLA_FV02_NRT00_C-20220718T224508Z.nc.gz'];
filename = ['https://thredds.aodn.org.au/thredds/dodsC/IMOS/OceanCurrent/GSLA/NRT00/2022/IMOS_OceanCurrent_HV_20220722T060000Z_GSLA_FV02_NRT00_C-20220724T224221Z.nc.gz'];
%ssh =load([science_path,'\DATA\SATELLITE\SSH\20220715.mat']);
disp(filename);

%---------------------------------------------------------------------------------------------------------------------------
% get SSH
%-------------
SSH.LONGITUDE = ncread(filename,'LONGITUDE');
SSH.LATITUDE = ncread(filename,'LATITUDE');
SSH.UCUR = ncread(filename,'UCUR');
SSH.VCUR = ncread(filename,'VCUR');
SSH.TIME = ncread(filename,'TIME')+datenum(1985,01,01);
[SSH.X, SSH.Y] = meshgrid(SSH.LONGITUDE,SSH.LATITUDE);
SSH.UCUR = SSH.UCUR';
SSH.VCUR = SSH.VCUR';

%%%%%%%%%5 OR FrOM BENOIT
% lon_ssh=ncread([science_path,'\DATA\Satellite\nrt_merged_alti_East_Australian_Current_2021_05.nc'],'longitude') ;
% lat_ssh=ncread([science_path,'\DATA\Satellite\nrt_merged_alti_East_Australian_Current_2021_05.nc'],'latitude') ;
% ugos_ssh=ncread([science_path,'\DATA\Satellite\nrt_merged_alti_East_Australian_Current_2021_05.nc'],'ugos') ;
% vgos_ssh=ncread([science_path,'\DATA\Satellite\nrt_merged_alti_East_Australian_Current_2021_05.nc'],'vgos') ;
% sla_ssh=ncread([science_path,'\DATA\Satellite\nrt_merged_alti_East_Australian_Current_2021_05.nc'],'sla') ;
% %sla_ssh=ncread('\\data.investigator.csiro.au\voyages\current\science\DATA\Satellite\nrt_merged_alti_East_Australian_Current_2021_05.nc','sla') ;
% 

%% PLOT SST
figure
set(gcf, 'position',[50   50   600   600],'PaperPositionMode','auto' , 'Color', 'w')
hold on
pcolor(sst.lns,sst.lts,sst.sst1); shading interp        % contourf(x,y,SST,50,'linestyle','none')
caxis([18 23])
%colormap(cmocean('thermal'))
%colormap(rainbow)
xl = [153 156];
yl = [-29.5 -25];
colormap(cmocean('thermal'))
hh=colorbar;  set(get(hh,'ylabel'),'String', 'Sea Surface Temperature (^oC)','Fontsize',10);
ylabel('Latitude (^oN)','Fontsize',10);
xlabel('Longitude (^oE)','Fontsize',10);
title(['SST: ' datestr(datevec(date_sst,'yyyymmddHH'),'dd mmm yyyy HH') 'h, SSH: ' datestr(SSH.TIME,'dd mmm yyyy HH') 'h'],'Fontsize',10);
xlim(xl);  ylim(yl);
contour(lon_bat,lat_bat,height_bat',[ -2000 -1000 -200],'color',[0.1 0.1 0.1]);
contour(lon_bat,lat_bat,height_bat',[-400 -500],'color',[0.8 0.8 0.8]);
plot(eaccoast(:,1),eaccoast(:,2),'-k','LineWidth',1);
% Ratio lat / lon
Y=get(gca,'ylim');
set(gca,'dataaspectratio',[1 cos((Y(2)+Y(1))/2/180*pi) 1])

% export_fig(['../FIGURES/fig_SST_' datestr(time_sst,'yyyymmdd') '_article.png'],'-m2', '-q101')
% export_fig(['../FIGURES/fig_SST_' datestr(time_sst,'yyyymmdd') '_article.pdf'],'-painters', '-r5', '-q101')

%u_SSH(find(lat_SSH==-33)+1,find(lon_SSH==150)+1)=1; v_SSH(find(lat_SSH==-33)+1,find(lon_SSH==150)+1)=0;
% quiver(lon_SSH(1:2:end),lat_SSH(1:2:end),u_SSH(1:2:end,1:2:end),v_SSH(1:2:end,1:2:end),1,'color','k')
% text(150,-33,' 1 m s^{-1}   ','HorizontalAlignment','left','Fontsize',10,'color','k')


hold on
quiver(SSH.LONGITUDE,SSH.LATITUDE,SSH.UCUR,SSH.VCUR,3,'k');
% OR
%quiver(lon_ssh,lat_ssh,ugos_ssh(:,:,end)',vgos_ssh(:,:,end)',3,'k');



%%%%%%%%%%%%%5 Add moorings
mooringlons = [153.90, 154.002, 154.134, 154.286, 154.639, 155.30230]
mooringlats = [-27.325, -27.316, -27.281, -27.244, -27.202, -27.1026]
plot(mooringlons,mooringlats,'o','MarkerFaceColor','k','MarkerSize',10)
% Triaxus
p = plot([154.61 153.72],[-27 -27.14],'LineWidth',2,'Color','k');
%p = plot([154.61 153.76],[-27.00 -27.12],'LineWidth',2,'Color','g');
%p = plot([155.3 153.8],[-27.42 -27.7],'LineWidth',1,'Color','g');
p = plot([155.3 154.45],[-27.42 -27.567],'LineWidth',2,'Color','k');
p = plot(153.72053881,-27.24306728,'LineWidth',2,'Color','r');

%19th of July: bad weather CTD + rapic cast line
% ctds_lon = [153.59776860, 153.63890122, 153.68022060, 153.72053881, 153.76078592, 153.80206560, 153.84237918, 153.88238850, 153.92300721,153.96270780,154.00195206]
% ctds_lat = [-27.26301150, -27.25593521, -27.24965017, -27.24306728, -27.23672612, -27.23068926, -27.22338391, -27.21756669, -27.21066719,-27.20457939,-27.19808575]
ctds_lon = [153.59776860, 153.63890122, 153.68022060, 153.72053881, 153.76078592]
ctds_lat = [-27.26301150, -27.25593521, -27.24965017, -27.24306728, -27.23672612]
plot(ctds_lon,ctds_lat,'o','MarkerFaceColor','r','MarkerSize',5)

% SBP line
p = plot([154.050859 154.068007],[-27.556572 -28.371804],'LineWidth',2,'color',[0.4 0.4 0.4]);
% Northern end of Canyon
canyonE_lons = [154.580035, 154.541667, 154.505856]
canyonE_lats = [-28.515384, -28.629365, -28.731207]
plot(canyonE_lons,canyonE_lats,'-og','MarkerFaceColor','g','MarkerSize',2)
canyonMid_lons = [154.247828, 154.216911, 154.184715]
canyonMid_lats = [-28.466660, -28.527136, -28.591823]
plot(canyonMid_lons,canyonMid_lats,'-og','MarkerFaceColor','g','MarkerSize',2)
canyonW_lons = [154.009276, 153.973671, 153.946479]
canyonW_lats = [-28.294091, -28.330575, -28.358670]
plot(canyonW_lons,canyonW_lats,'-og','MarkerFaceColor','g','MarkerSize',2)

canyonmouth_lons = [154.541667, 154.216911, 153.973671]
canyonmouth_lats = [-28.629365, -28.527136, -28.330575]
plot(canyonmouth_lons,canyonmouth_lats,'-og','MarkerFaceColor','g','MarkerSize',2)

plot(154.25,-26.94,'o','MarkerFaceColor','r','MarkerSize',5)
% p = plot([154.25 153.8525],[-26.94 -27.006],'LineWidth',2,'color',[0.4 0.4 0.4]);
p = plot([154.3462 153.8525],[-26.934 -27.006],'LineWidth',2,'color',[0.4 0.4 0.4]);
p = plot([154.3462 153.8146],[-26.944 -27.015],'LineWidth',2,'color',[0.4 0.4 0.4]);


%%% Save
export_fig([science_path,'FIGURES/fig_matlab_satellite_SST_mat_' date_sst '_triaxusline1_2_CTDsSBP_canyon.png'],'-m2', '-q101')
% export_fig([science_path,'FIGURES/fig_matlab_satellite_SST_mat_' date_sst '_triaxusline1.png'],'-m2', '-q101')
%export_fig([science_path,'FIGURES/fig_SST4h_Griffin.pdf'],'-painters', '-r5', '-q101')


