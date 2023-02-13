function [A,t,MassFluxes]=BaratariaArea_final(X,x,f,eustatic,M_f,M_sm,deepsubs,nyr,t,dt,f_rec,PostDamSedFlux,CouvillionTable1)
% sediment mass balance mode to relate sediment fluxes in Barataria Basin 
% to land area change. Edmonds et al. 

%Governing Equation:
%dA/dt = Anew - Aero - Arslr
%A_new = rate new marsh area creation
%A_ero = rate or marsh area erosion
%A_rslr= rate of marsh area loss by relative sea level rise

%% inputs
MinMass= X(1);   %mass of mineral sediment accumulated in Barataria in kg/km2/yr
OrgMass= X(2);   %mass of mineral sediment accumulated in Barataria in kg/km2/yr
SS = X(3);       %mean value of shallow subsidence over the whole delta in m/yr
d_avg = X(4);    %mean water depth in Barataria (m)
f_ero = X(5);    %fraction of landloss due to edge erosion (--)
InitArea = X(6); %Initial area of Barataria (km2)

%% constants
rhoWet=1150; %density (kg/m3) of sediment eroded from the wetland, from Bomer et al., (2019) as reported in Sanks et al., (2020)
rhoSurf=290; %density (kg/m3) for freshly deposited surficial sediment, from Sanks et al., (2020)

%% Pre-allocate variables
z_n=zeros(length(t)+1,1);
f_sub=zeros(length(t)+1,1);
A_ero=zeros(length(t),1);
M_or=zeros(length(t),1); 
M_rec=zeros(length(t),1);
M_oc=zeros(length(t),1);
M_e=zeros(length(t),1);
V_dep=zeros(length(t),1);
V_rslr=zeros(length(t),1);
sedBudget=zeros(length(t),1);
A_rslr=zeros(length(t),1);
A_new=zeros(length(t),1);
nRecEff=numel(f_rec);
A=zeros(length(t),nRecEff);
shallowsubs=ones(size(eustatic)).*SS;
%% calculate variables from input
TotMass=OrgMass+MinMass; %total organic + mineral depositional mass on Barataria, kg/km2/yr
OrgFrac=OrgMass/TotMass; %fraction of depositional mass that is organic

%Calculate f_Oc depending on value of f_rec. 
BaratariaArea_Sanks=3083; %area of barataria when Sanks et al., (2020) reported from their paper
AvgA=mean(CouvillionTable1.Area_km2(14:20)); %average Barataria area from 2006 to 2015 using Table 1 in Couvillion et al., (2017) (km2)
Aero_decade=(f_ero*AvgA)*1000^2; %Area of erosion over the decade of measurement by Sanks et al. (m2/yr) 
Rec_Mass=f_rec.*(Aero_decade*d_avg*rhoWet); %proportion of land loss that is recylced back onto the marsh surface for different efficiencies (kg/yr)
Oc_mass=(MinMass*BaratariaArea_Sanks)-Rec_Mass; %mineral mass accumulation in Barataria without recycled component (kg/yr). We consider this to be the offshore mass transported onshore -- so-called ocean mass
if Oc_mass < 0
    f_Oc = 0;
else
    f_Oc=(Oc_mass)./PostDamSedFlux; %Equation 7 from Edmonds et al. that gives the retention factor from oceans. 
end
%% Solve time-dependent change in area of Barataria Basin

A(1,:)=InitArea*1000^2; %Initial area (m2) of Barataria Basin
R=eustatic+shallowsubs+deepsubs;

%Compute the change in area for a given scenario for different recycling efficiency  
    for i = 1:length(t)-1     
        M_or(i,1)=TotMass*OrgFrac*(A(i,1)/1000^2); %Equation 3 from Edmonds et al.--mass of organic sediment added to marsh per year (kg/yr)
        A_ero(i,1)=A(i,1)*f_ero; %Equation 5 from Edmonds et al.--marsh edge erosion rate (m2/yr)  
        M_e(i,1) = A_ero(i,1)*d_avg*rhoWet; %rate of mass erosion of sediment from marsh edges (kg/yr)
        M_rec(i,1)=f_rec(1).*(M_e(i,1)); %Equation 4 from Edmonds et al.--rate of sediment mass erosion from marsh edges that is recylced back to marsh (kg/yr)
              
        M_oc(i,1)=(M_sm(i)-M_f(i)).*f_Oc(1); %Equation 6 from Edmonds et al., the oceanic contribution to deposition (kg/yr)
        V_dep(i,1)=(M_f(i)+M_rec(i,1)+M_or(i,1)+M_oc(i,1))*(1/rhoSurf); %Equation 8 from Edmonds et al.--volumetric rate (m3/yr) of contribution from fluvial overbank deposition (M_f), organic deposition (M_or), marsh edge recylcine (M_r), and ocean derived (M_oc) 
        
        V_rslr(i,1)=A(i,1)*(R(i)); %Equation 9 from Edmonds et al.--volume needed to keep up with relative sea level rise (m3/yr)

        sedBudget(i,1)=V_dep(i,1)-V_rslr(i,1); 
           
        if sedBudget(i,1)<0 %sed deficit
           z_n(i,1)=(R(i)-(V_dep(i,1)/A(i,1)))*dt + z_n(i-1,1); %Equation 11 from Edmonds et al.--relative sea level change minus sediment deposition yields increase in R relative to land surface. Add previous timestep to track total increase and land flooded       
           f_sub(i,1) = interp1(x,f,z_n(i,1),'nearest','extrap'); %get value of proportion of area flooded cdf of elevation (x and f) for that z_n      
           A_rslr(i,1)=A(i,1)*(f_sub(i,1)-f_sub(i-1,1)); %Equation 12 from Edmonds et al.--land submerged by RSLR (m2)
           A(i+1,1)=A(i,1)-A_rslr(i,1)-A_ero(i,1)*dt; %next marsh area: land loss by subsidence and erosion      
        elseif sedBudget(i,1)>=0 %sed surplus
           z_n(i,1)=0; %no change in relative sea level since sedimentation keeps up
           f_sub(i,1)=0;
           A_new(i,1)=sedBudget(i,1)*rhoSurf/rhoWet/d_avg; %Equation 10 from Edmonds et al.--determine areal growth (m2/yr)
           A(i+1,1)=A(i,1)+A_new(i,1)*dt-A_ero(i,1)*dt; %next marsh area
        end
        
    end
    MassFluxes=[M_or M_rec M_f M_oc];
   


