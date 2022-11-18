%%This code calculate the area change in Barataria Basin from 1890 to 2020
%%under 5 different scenarios defined in Edmonds et al. The code is
%%presently set up to run with one recylcing scenario (frec = 0.75). This
%%code reproduces Fig. 4A. To reproduce all of figure 4, save the output
%%from this code and run the script Figure4.m. All references numbers refer
%%to the main text

clear all
close all
warning off

%% setup model constants
startYear = 1890; %beginning of simulation
endYear = 2020; %end
nyr = endYear - startYear; %number of years simulated
dt = 1; % time step in years
t = ((1:dt:nyr) - 1) + startYear; % total years simulated
recEff = 0.75; %0:0.25:1; %recycling efficiency scenarios (f_rec)
nrecEff = numel(recEff);

DS = 0.0032; %deep subsidence value (m/yr) from ref. 14 in Edmonds et al.
PredamSedFlux = 348 * 1000000000; %pre-dam sediment flux of the Mississippi river, converted from megatons/yr to kg/yr (from ref. 20 in Edmonds et al.) 
PostdamSedFlux = 91 * 1000000000; %post-dam sediment flux of the Miss River, converted from megatons/yr to kg/yr (from ref. 13 in Edmonds et al.)
kg_overbank = 32 * 1000000000; %this is in kg of sediment deposited overbank from D3D simulation, converted from megatons/yr to kg/yr. calculated with CalculateOverbankDeposition.m
YrDamBuild = 1946 - startYear; %1946 is the midpoint year of most dam building from ref. 20 in Edmonds et al.
YrLeveeBuild = 1935 - startYear; %1935 is the midpoint year of most levee building 
postdam_reduction = 0.66; %this is the proportional reduction in sediment deposition following dam construction from D3D simulations
natural_wave_erosion = 0.30; %fraction of all erosion from 1932 to 1990 attributed to 'natural' edge erosion from wave (from ref. 28 in Edmonds et al.) 

%% Load land loss and elevation data
load('CouvillionTable1.mat'); %from ref. 5 in Edmonds et al. main text
areaBar = CouvillionTable1.Area_km2;
areaBar_error = CouvillionTable1.SE_km2; 
yr = CouvillionTable1.Year;
load('CouvillionTable3.mat')
load('Barataria_elevations_gt0_FINAL', 'f_gt0', 'x_gt0') %data from ref. 44 in Edmonds et al. 
x = x_gt0(2:end); %ecdf duplicates the first point, need to remove for interp1
f = f_gt0(2:end);
load('Barataria_elevations_lt0_FINAL', 'lt0') %data from ref. 44 in Edmonds et al. 
z_lt0 = [lt0.grid_code];

%% Load global mean sea level (GMSL)
table = readtable('global_basin_timeseries.xlsx', 'Sheet', 'Subtropical North Atlantic'); %data from ref. 22 in Edmonds et al.
sl = table.ObservedBasin_meanSeaLevel_mean_; %average sea level through time for North Atlantic (mm)
sl_time = table.Var1; %year of sea level measurement
p = polyfit(sl_time, sl, 3); %3rd-order polynominal fit to the data
sl_fit = polyval(p, startYear:endYear)'; %evaluate that polynomial over the simulation time
sl_fit_zero1890 = sl_fit - sl_fit(1); %set sea level to zero in 1890
sl_rise = diff(sl_fit_zero1890); %yearly rate of sea level rise (mm/yr)
eustatic = sl_rise / 1000; %convert to m/yr

%% Load Kolker et al. 2011 data (ref. 11)
load Kolker2011_polyfit %from ref. 11 in Edmonds et al.
intpYr = 1949:2003; %years for interpolation
intpsubs = feval(fittedmodel, intpYr)'; %evaulate interpolated data at points specified by intpYr
intpsubs = intpsubs / 1000; %convert mm/yr to m/yr

%% initial guess of Monte Carlo Variables
MinMass_mean = (6.92 / 3083) * 1000000000; %mass of mineral sediment accumulated in Barataria in kg/km2/yr from ref. 13 in Edmonds et al. 
MinMass_SE = (3.56 / 3083) * 1000000000; %one standard error on TotSedMassBarataria in Barataria kg/km2/yr from ref. 13 in Edmonds et al.

OrgMass_mean = (2.91 / 3083) * 1000000000; %mass of mineral sediment accumulated in Barataria in kg/km2/yr from ref. 13 in Edmonds et al.
OrgMass_SE = (1.22 / 3083) * 1000000000; %one standard error on OrgMass in Barataria kg/km2/yr from ref. 13 in Edmonds et al.

SS_mean = 0.0068; %mean value of shallow subsidence over the whole delta in m/yr from Table 1 in ref. 14 in Edmonds et al.
SS_SE = 0.0079 / (sqrt(274)); %one standard error in m/yr from Table 1 in ref. 14 in Edmonds et al.--Standard Error is calculated as standard deviation divided by sqrt of num observations

d_mean = abs(mean(z_lt0(z_lt0 < -0.2))); %mean water depth in Barataria for all depths less than -0.2. This is done to avoid averaging over the shallow intertidal platform depths since we want marsh height/thickness
d_std = abs(std(z_lt0(z_lt0 < -0.2))); %standard deviation
N = length(z_lt0(z_lt0 < -0.2));
d_SE = d_std / sqrt(N)*300; %standard error of the mean, use 300 to increase error from the otherwise absurbly small value

[f_ero_mean, f_ero_SE] = f_ero_mean_error(CouvillionTable1, CouvillionTable3,natural_wave_erosion); %calculate fraction of landloss due to edge erosion with error propagation

Aknot = 3200; %guess of the initial Barataria Area, km2
Aknot_SE = 2*areaBar_error(1); %increase standard error to 95% conf. int with 1.96 multiplier and assume standard error from the first time step measured is good approx for 1890 

%% Load initial Guesses for Monte Carlo
X_mean = [MinMass_mean, OrgMass_mean, SS_mean, d_mean, f_ero_mean, Aknot];
SE = [MinMass_SE, OrgMass_SE, SS_SE, d_SE, f_ero_SE, Aknot_SE];

%% number of Monte Carlo samples
NumSamples = 5e2;
cnt = 0;
%% MAIN LOOP--calculate five different model scenarios doing 1000 monte carlo simulations of each
for k = 1:5
    
    M_sm = zeros(length(t), 1); %Sediment flux of the Mississippi River
    M_f = zeros(length(t), 1); %overbank flux of sediment deposited in Barataria
    MassFluxesAvg=zeros(length(t),4);

    if k == 1 %pre-dam scenario
        M_f(1:end) = kg_overbank; %The depositional flux (kg/yr) of overbank sediment flux from the Mississippi River to Barataria Basin from Delft3D
        M_sm(1:end) = (PredamSedFlux); %Predam sediment flux in the Mississippi River (kg/yr)
        deepsubs = ones(size(eustatic)) .* DS;
    end
    if k == 3 %dams-only scenario
        M_f(1:YrDamBuild/dt) = kg_overbank; %The depositional flux (kg/yr) of overbank sediment flux from the Mississippi River to Barataria Basin from Delft3D
        M_f(YrDamBuild/dt:length(t)) = kg_overbank * postdam_reduction; %effect of dams is to reduce M_f
        M_sm(1:YrDamBuild/dt) = (PredamSedFlux ); %Predam sediment flux in the Mississippi River (kg/yr) UNITS?
        M_sm(YrDamBuild/dt:end) = (PostdamSedFlux); %Postdam sediment flux in the Mississippi River (kg/yr)
        deepsubs = ones(size(eustatic)) .* DS;
    end
    if k == 2 %levees-only scenario
        M_f(1:YrLeveeBuild/dt) = kg_overbank; %The depositional flux (kg/yr) of overbank sediment flux from the Mississippi River to Barataria Basin from Delft3D
        M_f(YrLeveeBuild/dt:length(t)) = 0; %after levees are built M_f is 0
        M_sm(1:end) = (PredamSedFlux); %Predam sediment flux in the Mississippi River (kg/yr)
        deepsubs = ones(size(eustatic)) .* DS;
    end
    if k == 4 %dams and levee scenario
        M_f(1:YrLeveeBuild/dt) = kg_overbank; %The depositional flux (kg/yr) of overbank sediment flux from the Mississippi River to Barataria Basin from Delft3D
        M_f(YrLeveeBuild/dt:length(t)) = 0; %after levees are built M_f is 0
        M_sm(1:YrDamBuild/dt) = (PredamSedFlux); %Predam sediment flux in the Mississippi River (kg/yr) UNITS?
        M_sm(YrDamBuild/dt:end) = (PostdamSedFlux); %Postdam sediment flux in the Mississippi River (kg/yr)
        deepsubs = ones(size(eustatic)) .* DS;
    end
    if k == 5 %resource extraction scenario
        M_f(1:YrLeveeBuild/dt) = kg_overbank; %The depositional flux (kg/yr) of overbank sediment flux from the Mississippi River to Barataria Basin from Delft3D
        M_f(YrLeveeBuild/dt:length(t)) = 0;
        M_sm(1:YrDamBuild/dt) = (PredamSedFlux); %Predam sediment flux in the Mississippi River (kg/yr) UNITS?
        M_sm(YrDamBuild/dt:end) = (PostdamSedFlux); %Postdam sediment flux in the Mississippi River (kg/yr)

        % resource extraction affects DS, sub in DS values from ref. 11
        difference = (intpsubs - DS)';
        difference(difference < 0) = 0;
        index = intpYr - startYear;
        deepsubs = ones(size(eustatic)) .* DS;
        deepsubs(index) = deepsubs(index) + (difference);
    end

     for j = 1:nrecEff %number of recylcing scenarios
        f_rec = recEff(j); %recycling fraction
        for n = 1:NumSamples
            accept_sample=0;
            X = ones(size(X_mean)) .* -1; %setup dummy variable to track if monte carlo selected values are positive or negative. Initialize as all negative.

            while sum(X < 0) > 0 %run this while loop until there are five X values that are all positive
                ind = X < 0; %index locations where Xnew < 0
                r = normrnd(X_mean(ind), SE(ind)); %select a random value for mean given the standard error
                X(ind) = r; %reset Xnew choice to previous choice X
            end          
           
                [A, t,MassFluxes] = BaratariaArea_final(X, x, f, eustatic, M_f, M_sm, deepsubs, nyr, t, dt, f_rec,PostdamSedFlux,CouvillionTable1); %calculate Barataria area by solving equation 1 in Edmonds et al.
                MassFluxesAvg=MassFluxes+MassFluxesAvg;
                A = A ./ 1000^2;

            cnt = cnt + 1;
            M(:, cnt) = X'; %set of randomly selected values for this prediction of A
            Ahats(:, cnt) = A; %save prediction for A through time in matrix Ahats

            disp(['Monte Carlo iteration #', num2str(n), ' of ', num2str(NumSamples)])
        end

        linecolors = [0, 0.6, 0.77; 0.89, 0, 0.23; 0.4940, 0.1840, 0.5560; 0.9290, 0.6940, 0.1250; 0.8500, 0.3250, 0.0980];
        Ahats_mean = mean(Ahats, 2);
        Ahats_SE = std(Ahats') ./ sqrt(cnt);
        CI95 = tinv([0.025, 0.975], cnt-1);
        yCI95 = bsxfun(@times, Ahats_SE, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of hxa
        figure(1)
        hold on
        plot(t, Ahats_mean, 'Color', linecolors(k, :), 'LineWidth', 2)
        %plot(t, Ahats, 'Color', linecolors(k, :), 'LineWidth', 0.25)
        patch([t, fliplr(t)], [Ahats_mean' + yCI95(1, :), fliplr(Ahats_mean'+yCI95(2, :))], linecolors(k, :), 'FaceAlpha', 0.5, 'EdgeColor', 'none')
        cnt = 0;
        if k == 1
            A_predam = Ahats;
            A_predam_SE = Ahats_SE;
            M_predam = M;
            MassFluxesAvg_predam=MassFluxesAvg./n;
        end
        if k == 2
            A_dams = Ahats;
            A_dams_SE = Ahats_SE;
            M_dams = M;
            MassFluxesAvg_dam=MassFluxesAvg./n;
        end
        if k == 3
            A_levees = Ahats;
            A_levees_SE = Ahats_SE;
            M_levees = M;
            MassFluxesAvg_levees=MassFluxesAvg./n;
        end
        if k == 4
            A_damslevees = Ahats;
            A_damslevees_SE = Ahats_SE;
            M_damslevees = M;
            MassFluxesAvg_damslevees=MassFluxesAvg./n;
        end
        if k == 5
            A_damsleveesresource = Ahats;
            A_damsleveesresource_SE = Ahats_SE;
            M_resource = M;
            MassFluxesAvg_damsleveesresource=MassFluxesAvg./n;
        end
        clear Ahats 
        MassFluxesAvg = 0;

    end
end
errorbar(yr,areaBar,2.*areaBar_error,'ko','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',6,'LineWidth',0.5)

ax = gca;      
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlim([startYear endYear])
ylim([2500 5000])
yticks([2500 3000 3500 4000 4500 5000 )
ax.YAxis.MinorTick = 'on';
ax.XAxis.MinorTick = 'on';
ylabel('Barataria Area (km^2)')
xlabel('Calendar Year')
box on
set(gca,'linewidth',1)
ax.TickLength = [.02, .02];