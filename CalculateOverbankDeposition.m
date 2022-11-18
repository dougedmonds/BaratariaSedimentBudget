clear all
close all

%This code uses Delft3D FM results to recreate figure 3 from Edmonds et al.
%and also calculate the yearly deposition for pre and post dam conditions.

fName='FlowFM_map'; %d3d results unmodified
areaName='mesh2d_flowelem_ba'; %area of each grid
diffName='mesh2d_dpsed';
massName='mesh2d_bodsed';
load('MRD_allexport','allexport'); %latitude longitude of all grid cell centres, as exported from arcmap
load('MRD_flpl_lalo','clipFloodplain'); %floodplain points exported from arcmap
iFp=ismember(allexport,clipFloodplain,'rows'); %find indices of rows that contain floodplain coordinate

for j =1:2 %loop through predam simulations, then postdam
    if j == 1
        cd predam_simulations
        dirnames=dir('*Q.');
    else
        cd postdam_simulations
        dirnames=dir('*Q.');
    end

for i = 1:length(dirnames)
         name=dirnames(i).name;
         cd([name])

         cd Project1.dsproj_data
         cd FlowFM
         cd output

        %get area of each gridcell
        ncName=strcat(fName,'.nc'); %name of each run
        gridAreaAll=ncread(ncName,areaName); %extract initial topo map
        gridAreaFp=gridAreaAll;
        gridAreaFp(iFp==0)=0; %gridcell area, zeros where it's not floodplain
        totalAreaFp=sum(gridAreaFp); %total floodplain area

        %extract mass deposition
        temp=ncread(ncName,massName);
        for sedfrac=1:4
            dMass(:,sedfrac)=squeeze(temp(sedfrac,:,end)-temp(sedfrac,:,1));
            dMass(iFp==0,sedfrac)=0; %exclude everything that's not floodplain
        end
        
        massLocal=dMass.*gridAreaFp; %mass gained in each grid cell (m3)
        massFrac=sum(massLocal); %sum mass deposited in all grid cells
        massRun(:,i)=sum(massFrac)/2;
        
        cd ..
        cd .. 
        cd ..
        cd ..   
        
end

cd ..

run=20:5:40;
%depositional volume
figure(1)
hold on
plot(run*1000,massRun/1000000000,'-o')
xlabel('Peak flood discharge (x10^3 m^3/s)')
ylabel('Deposited mass (MT)')
hold off

%%
%convert individual floods to a mean aggradation rate per year from USGS
%gauge 07374000
load('floodRecurrence_data.mat');
pFl(pFl<20000)=[];

%Slope between each discharge via linear interpolation
d1=(massRun(2)-massRun(1))/5000;
d2=(massRun(3)-massRun(2))/5000;
d3=(massRun(4)-massRun(3))/5000;
d4=(massRun(5)-massRun(4))/5000;
d5=d4; %temporary, until 45k run works

nFloodMeasured=numel(pFl);
floodDep=zeros(nFloodMeasured,1);

for i=1:nFloodMeasured
    peakFlow=pFl(i);
    if peakFlow>=20000 && peakFlow<=25000
        floodDep(i)=massRun(1)+d1*(peakFlow-20000);
    end
    if peakFlow>25000 && peakFlow<=30000
        floodDep(i)=massRun(2)+d2*(peakFlow-25000);
    end
    if peakFlow>30000 && peakFlow<=35000
        floodDep(i)=massRun(3)+d3*(peakFlow-30000);
    end
    if peakFlow>35000 && peakFlow<=40000
        floodDep(i)=massRun(4)+d4*(peakFlow-35000);
    end
    if peakFlow>40000 && peakFlow<=45000
        floodDep(i)=massRun(5)+d5*(peakFlow-40000);
    end
    if peakFlow>45000
        fprintf('\nThis flood is larger than expected. (>45k)\n')
    end
%     hold on
%     plot(peakFlow,floodDep(i),'ro')
end

depYearly(j)=mean(floodDep)/1000000000; %yearly averaged depositional mass (MT/yr)

end





