%% Model path
c2vsim_path = fullfile('..','C2VsimV101','AreaWeightedDiv');
%% Read Groundwater Recharge
c2vsimTime = dateshift(datetime('20-Sep-1973'),'end','month',1:504)';
GBinfo = h5info(fullfile(c2vsim_path, "Results",'C2VSimFG_GW_ZBudget.hdf'));
colIDnames = h5read(GBinfo.Filename,...
    [GBinfo.Groups(1).Name GBinfo.Name GBinfo.Groups(1).Datasets(5).Name]);
colIDs = h5read(GBinfo.Filename,...
    [GBinfo.Groups(1).Name GBinfo.Name GBinfo.Groups(1).Datasets(6).Name]);
%%
DeepPerc = h5read(GBinfo.Filename,...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(5).Name]);
%%
DivLoss = zeros(size(DeepPerc));
tmp = h5read(GBinfo.Filename,...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(7).Name]);
divElIds = find(colIDs(:,17) ~= 0);
DivLoss(divElIds,:) = tmp;
%%
ByPassLoss = zeros(size(DeepPerc));
tmp = h5read(GBinfo.Filename,...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(1).Name]);
bypassElIds = find(colIDs(:,19) ~= 0);
ByPassLoss(bypassElIds,:) = tmp;
%%
c2vsimRch = DeepPerc + DivLoss + ByPassLoss;
%% Find the indices of the start and and averaging periods
idx_2000 = find(c2vsimTime == datetime('31-Mar-2000')):find(c2vsimTime == datetime('31-May-2000'));
idx_2015 = find(c2vsimTime == datetime('31-Mar-2015')):find(c2vsimTime == datetime('31-May-2015'));
%% Average recharge
Rch_2000 = (0.3048^3)*sum(c2vsimRch(:,idx_2000),2)./sum(c2vsimTime(idx_2000).Day); %m^3/day
Rch_2015 = (0.3048^3)*sum(c2vsimRch(:,idx_2015),2)./sum(c2vsimTime(idx_2015).Day); %m^3/day
%% Read element area
%ElemArea = h5read(GBinfo.Filename,...
%    [GBinfo.Groups(1).Name GBinfo.Name GBinfo.Groups(1).Datasets(12).Name]);
%ElemArea = ElemArea*(0.3048^2);
c2vsim_mesh = shaperead(fullfile(c2vsim_path,'..','gis_data','C2VSim_Elements_3310'));
ElemArea = zeros(length(c2vsim_mesh),1);
bc_elem = zeros(length(c2vsim_mesh),2);
for ii = 1:length(c2vsim_mesh)
    ElemArea(ii,1) = polyarea(c2vsim_mesh(ii,1).X(1:end-1), c2vsim_mesh(ii,1).Y(1:end-1));
    bc_elem(ii,:) = [mean(c2vsim_mesh(ii,1).X(1:end-2)) mean(c2vsim_mesh(ii,1).Y(1:end-2))];
end
%% Calculate recharge rates
Rch_2000 = 1000*365*Rch_2000./ElemArea;
Rch_2015 = 1000*365*Rch_2015./ElemArea;
%% 
figure()
clf
histogram(Rch_2000, 100, 'Normalization','probability','EdgeColor','none','DisplayName','Spring 2000')
hold on
histogram(Rch_2015, 100, 'Normalization','probability','EdgeColor','none','DisplayName','Spring 2015')
xlim([0 2000])
legend('Location','northeast')
grid on
xlabel('Groundwater Recharge [mm/year]')
ylabel('%')
% print -dpng -r300 RchHist
%%
gwl_data1 = readtable(fullfile('..','..','Box','cv-unsat','gwl_file_part_1.xlsx'));
gwl_data2 = readtable(fullfile('..','..','Box','cv-unsat','gwl_file_part_2.xlsx'));
gwl_data3 = readtable(fullfile('..','..','Box','cv-unsat','gwl_file_part_3.xlsx'));
%%
columns_to_keep = ["Var2","Var4","Var5","Var6","Var7","Var8"];
gwl_data = [gwl_data1(:,columns_to_keep) 
    gwl_data2(:,columns_to_keep)
    gwl_data3(:,columns_to_keep)];
gwl_data.Properties.VariableNames = {'Section', 'Date','Var5','Var6','Var7','Var8'};
gwl_data.Section = categorical(gwl_data.Section);
%%
gwl_data = gwl_data((gwl_data.Date >= datetime(2000,2,1) & gwl_data.Date <= datetime(2000,5,31)) | ...
                    (gwl_data.Date >= datetime(2015,2,1) & gwl_data.Date <= datetime(2015,5,31)),:);
%%
gwl_data.DGW = gwl_data.Var8 - (gwl_data.Var5 - gwl_data.Var6);
gwl_data(:,["Var5","Var6","Var7","Var8"]) = [];
gwl_data(isnan(gwl_data.DGW),:) = [];
%%
gst = readtable(fullfile('..','..','Box','cv-unsat','gst_file.xlsx'));
gst.SITE_CODE = categorical(gst.SITE_CODE);
%%
[Lia, Locb] = ismember(gwl_data.Section, gst.SITE_CODE);
gwl_data.Lat(Lia) = gst.LATITUDE(Locb(Lia));
gwl_data.Lon(Lia) = gst.LONGITUDE(Locb(Lia));
gwl_data = gwl_data(Lia,:);% Remove the records with zero coordinates
%%
trs_unique = unique(gwl_data.Section);
GWL = table(trs_unique,'VariableNames', {'Section'});
%%
for ii = 1:size(GWL,1)
    ind = find(gwl_data.Section == GWL.Section(ii));
    if ~isempty(ind)
        GWL.Lat(ii) = gwl_data.Lat(ind(1));
        GWL.Lon(ii) = gwl_data.Lon(ind(1));
        % find records for spring 2000
        iyr = year(gwl_data.Date(ind)) == 2000;
        GWL.DGW_2000(ii) = mean(gwl_data.DGW(ind(iyr)));
        % find records for spring 2015
        iyr = year(gwl_data.Date(ind)) == 2015;
        GWL.DGW_2015(ii) = mean(gwl_data.DGW(ind(iyr)));
    end
end
%%
CV_outline = shaperead(fullfile('..','C2VsimV101','gis_data','C2VSim_Outline_3310'));
%%
figure()
clf
plot(CV_outline.X, CV_outline.Y)
[xx,yy] = projfwd(projcrs(3310),GWL.Lat, -GWL.Lon);
hold on
plot(xx,yy,'.')
title('All records')
hold off
% print -dpng -r300 AllRecordsMap
%%
CV_outline_shape = polyshape(CV_outline.X, CV_outline.Y);
in_cv = CV_outline_shape.isinterior(xx,yy);
GWL(~in_cv,:) = [];
%%
figure()
clf
subplot(1,2,1);
plot(-GWL.Lon(~isnan(GWL.DGW_2000)), GWL.Lat(~isnan(GWL.DGW_2000)),'.')
title({'Records with Spring',['2000 DGW (' num2str(sum(~isnan(GWL.DGW_2000))) ')']})
axis equal
axis off
subplot(1,2,2);
plot(-GWL.Lon(~isnan(GWL.DGW_2015)), GWL.Lat(~isnan(GWL.DGW_2015)),'.')
title({'Records with Spring',['2015 DGW (' num2str(sum(~isnan(GWL.DGW_2015))) ')']})
axis equal
axis off
% print -dpng -r300 GWLSpringMaps
%%
C2VsimHead = readIWFM_headalloutput(fullfile(c2vsim_path,'Results','C2VSimFG_GW_HeadAll.out'), 30179, 4, 505, 1);
%%
sim_wtbl_2000 = 0.3048 * (C2VsimHead{319,2}(:,1) + C2VsimHead{320,2}(:,1) + C2VsimHead{321,2}(:,1))/3;
sim_wtbl_2015 = 0.3048 * (C2VsimHead{499,2}(:,1) + C2VsimHead{500,2}(:,1) + C2VsimHead{501,2}(:,1))/3;
%%
cv_nodes = readIWFM_Nodes(fullfile(c2vsim_path, 'Preprocessor','C2VSimFG_Nodes.dat'));
cv_gse = readIWFM_Stratigraphy(fullfile(c2vsim_path,'Preprocessor','C2VSimFG_Stratigraphy.dat'),30179, 4, 105);
cv_gse = 0.3048 * cv_gse(:,2);
%%
sim_dgw_2000 = cv_gse - sim_wtbl_2000;
sim_dgw_2015 = cv_gse - sim_wtbl_2015;
%%
[GWL.X_3310, GWL.Y_3310] = projfwd(projcrs(3310),GWL.Lat, -GWL.Lon); 
GWL.DGW_2000 = GWL.DGW_2000*0.3048;
GWL.DGW_2015 = GWL.DGW_2015*0.3048;

[lat,lon] = projinv(projcrs(26910), cv_nodes(:,1), cv_nodes(:,2));
[simX3310, simY3310] = projfwd(projcrs(3310),lat, lon);
%%
Fmeas2000 = scatteredInterpolant(GWL.X_3310(~isnan(GWL.DGW_2000)), ...
    GWL.Y_3310(~isnan(GWL.DGW_2000)), GWL.DGW_2000(~isnan(GWL.DGW_2000)), 'linear', 'nearest');
Fmeas2015 = scatteredInterpolant(GWL.X_3310(~isnan(GWL.DGW_2015)), ...
    GWL.Y_3310(~isnan(GWL.DGW_2015)), GWL.DGW_2015(~isnan(GWL.DGW_2015)), 'linear', 'nearest');

Fsim2000 = scatteredInterpolant(simX3310, simY3310, sim_dgw_2000, 'linear', 'nearest');
Fsim2015 = scatteredInterpolant(simX3310, simY3310, sim_dgw_2015, 'linear', 'nearest');
%%
DGW2000sim = Fsim2000(Fmeas2000.Points(:,1), Fmeas2000.Points(:,2));
DGW2015sim = Fsim2015(Fmeas2015.Points(:,1), Fmeas2015.Points(:,2));
%%
Fmeas2000sim = scatteredInterpolant(Fmeas2000.Points(:,1), Fmeas2000.Points(:,2), DGW2000sim, 'linear','nearest');
Fmeas2015sim = scatteredInterpolant(Fmeas2015.Points(:,1), Fmeas2015.Points(:,2), DGW2015sim, 'linear','nearest');
%%
SimMeas2000 = Fmeas2000sim(simX3310, simY3310);
SimMeas2015 = Fmeas2015sim(simX3310, simY3310);
%%
measError2000 = Fsim2000.Values - SimMeas2000;
measError2015 = Fsim2015.Values - SimMeas2015;
%%
MeasInterp2000 = Fmeas2000(simX3310, simY3310);
MeasInterp2015 = Fmeas2015(simX3310, simY3310);
%%
DGW2000 = MeasInterp2000 + measError2000;
DGW2015 = MeasInterp2015 + measError2015;