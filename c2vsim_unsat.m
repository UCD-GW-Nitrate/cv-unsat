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


