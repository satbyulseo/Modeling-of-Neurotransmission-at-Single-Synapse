%%first.m, for initializing;
%synapse_size = 600; % 200, 400, 600
%cleft_geo='c';
%pos=5; %5 for center release and 9 for edge release
%R='i';
Dglut2=0.4;% Dglut value for outer region of cleft? \n - zone2 
Dglut3=0.4;% Dglut value for inner region of cleft? \n - zone1
%Dglut2=0.1;%  High affinity edge
%Dglut3=0.1;% High affinity center
tin=0;
tf=10000;
%Ndt=100;
%av=10;
%inc=2e-4; %ms
Ndt=50000;
av=10;
inc=2e-3; %ms





%% Parameters for diffusion of glutamate
%%

%% Size of Synapse
synapse_size = input('Size of synapse, 200, 400, or 600 nm? \n');

if synapse_size == 600
    i = 60; j = 101; k = 101;
    cleft.i = 29:31; cleft.j = 21:81; cleft.k = 21:81;
    zone.i = 29:31; zone.j = 34:67; zone.k = 34:67;
    %Receptor Locations
    R_Loc = zeros(25,2);
    R_Loc(1,:)=[26,76];R_Loc(2,:)=[38,76];R_Loc(3,:)=[51,76];R_Loc(4,:)=[64,76]; R_Loc(5,:)=[76,76];
    R_Loc(6,:)=[26,64];R_Loc(7,:)=[38,64];R_Loc(8,:)=[51,64];R_Loc(9,:)=[64,64]; R_Loc(10,:)=[76,64];
    R_Loc(11,:)=[26,51];R_Loc(12,:)=[38,51];R_Loc(13,:)=[51,51];R_Loc(14,:)=[64,51];R_Loc(15,:)=[76,51];
    R_Loc(16,:)=[26,38];R_Loc(17,:)=[38,38];R_Loc(18,:)=[51,38];R_Loc(19,:)=[64,38];R_Loc(20,:)=[76,38];
    R_Loc(21,:)=[26,26];R_Loc(22,:)=[38,26];R_Loc(23,:)=[51,26];R_Loc(24,:)=[64,26];R_Loc(25,:)=[76,26];
    %R_Loc(26,:)=[51,44];R_Loc(27,:)=[51,32];R_Loc(28,:)=[57,44];R_Loc(29,:)=[70,32];R_Loc(30,:)=[51,21];
    %R_Loc(32,:)=[73,29];R_Loc(34,:)=[60,41];R_Loc(35,:)=[51,41];
    %add the edge of the end for straight line

elseif synapse_size == 400
    i = 60; j = 61; k = 81;
    cleft.i = 29:31; cleft.j = 21:61; cleft.k = 21:61;
    zone.i = 29:31; zone.j = 30:52; zone.k = 30:52 ;
 
    %Receptor Locations
   
    R_Loc = zeros(25,2);
    R_Loc(1,:)=[25,57];R_Loc(2,:)=[33,57];R_Loc(3,:)=[41,57];R_Loc(4,:)=[49,57]; R_Loc(5,:)=[57,57];
    R_Loc(6,:)=[25,49];R_Loc(7,:)=[33,49];R_Loc(8,:)=[41,49];R_Loc(9,:)=[49,49]; R_Loc(10,:)=[57,49];
    R_Loc(11,:)=[25,41];R_Loc(12,:)=[33,41];R_Loc(13,:)=[41,41];R_Loc(14,:)=[49,41];R_Loc(15,:)=[57,41];
    R_Loc(16,:)=[25,33];R_Loc(17,:)=[33,33];R_Loc(18,:)=[41,33];R_Loc(19,:)=[49,33];R_Loc(20,:)=[57,33];
    R_Loc(21,:)=[25,25];R_Loc(22,:)=[33,25];R_Loc(23,:)=[41,25];R_Loc(24,:)=[49,25];R_Loc(25,:)=[57,25];
    %R_Loc(26,:)=[41,37];R_Loc(27,:)=[41,29];R_Loc(28,:)=[45,37];R_Loc(29,:)=[53,29];R_Loc(30,:)=[41,21];%add the edge of the end for straight line
    %R_Loc(32,:)=[55,27];R_Loc(33,:)=[45,25];R_Loc(34,:)=[47,35];R_Loc(35,:)=[41,35];
    
elseif synapse_size == 200
    i = 60; j = 61; k = 61;
    cleft.i = 29:31; cleft.j = 21:41; cleft.k = 21:41;
    zone.i = 29:31; zone.j = 26:36; zone.k = 26:36;
    
    %Receptor Locations
    R_Loc = zeros(32,2);
    R_Loc(1,:)=[23,39];R_Loc(2,:)=[27,39];R_Loc(3,:)=[31,39];R_Loc(4,:)=[35,39]; R_Loc(5,:)=[39,39];
    R_Loc(6,:)=[23,35];R_Loc(7,:)=[27,35];R_Loc(8,:)=[31,35];R_Loc(9,:)=[35,35]; R_Loc(10,:)=[39,35];
    R_Loc(11,:)=[23,31];R_Loc(12,:)=[27,31];R_Loc(13,:)=[31,31];R_Loc(14,:)=[35,31];R_Loc(15,:)=[39,31];
    R_Loc(16,:)=[23,27];R_Loc(17,:)=[27,27];R_Loc(18,:)=[31,27];R_Loc(19,:)=[35,27];R_Loc(20,:)=[39,27];
    R_Loc(21,:)=[23,23];R_Loc(22,:)=[27,23];R_Loc(23,:)=[31,23];R_Loc(24,:)=[35,23];R_Loc(25,:)=[39,23];
    %R_Loc(26,:)=[31,29];R_Loc(27,:)=[31,25];R_Loc(28,:)=[33,29];R_Loc(29,:)=[37,25];R_Loc(30,:)=[31,21];R_Loc(31,:)=[21,23];
    %R_Loc(32,:)=[38,24];R_Loc(33,:)=[33,23];R_Loc(34,:)=[34,28];R_Loc(35,:)=[31,28];
else
    disp('Size must either be 200, 400, or 600 nm. \n');
end
%%

%% Geometry of cleft, reduced edge or middle?
cleft_geo = input('Reduce the cleft on outer edge (e), center (c), or flat (f)?\n','s');
%%

%% Release location
pos = input('Release over which receptor, numeric value only. \n'); 
if (pos == 13 && cleft_geo == 'c') || (pos == 25 && cleft_geo == 'e')
    Release_loc(1) = cleft.i(1)+1; 
    Release_loc(2) = R_Loc(pos,1); 
    Release_loc(3) = R_Loc(pos,2);
else
    Release_loc(1) = cleft.i(1); 
    Release_loc(2) = R_Loc(pos,1); 
    Release_loc(3) = R_Loc(pos,2);
end
            
        
%%


%%



%% Diffusion coefficients
%Set Diffusion constants, Dglut1 is outside the cleft and Dglut2 is inside
%the cleft. dt = 0.2 microseconds. 
Dglut1=0.75;    %Outside the cleft
%Dglut2 = input('Dglut value for outer region of cleft? \n');
%Dglut2=0.40;    %Outer region of cleft
%Dglut3 = input('Dglut value for inner region of cleft? \n');
%Dglut3=0.2;    %Innermsot region of cleft

%dx=0.16;
%dx=0.08;
%dx=0.04; 
%dx=0.02;
dx=0.01; %0.01 micrometer or 10 nm 
%dx=0.005; 
%dx=0.0025; 

%dt=2e-5; %2e-5 ms or 0.02 microseconds

%dt = 0.01;
%dt = 0.005;
%dt = 0.0025;
%dt = 0.00125;
%dt = 0.000625;
%dt = 0.0003125;
%dt = 0.00015625;
%dt = 0.000078125;
%dt = 5.8594e-05;

%start stable
%dt = 3.90625E-05;
%dt = 1.95313E-05;
%dt = 4e-5;
%dt=2e-5;
%dt=1e-5;
%dt = 5e-6;
dt=2e-5; 


D1=(Dglut1*dt)/(dx)^2
D2=(Dglut2*dt)/(dx)^2
D3=(Dglut3*dt)/(dx)^2
Ddiff12=Dglut1/Dglut2;
Ddiff23=Dglut2/Dglut3;
%%

%Glutamate diffusion with multiple values within the synaptic cleft. 
%Cleft Glutamate Concentration Simulations

OK = 1;
global S C

R = input('Instantaneous, I, or Vesicle release, V?\n','s');
switch lower(R)
    case 'i'
        initial = 4000;
    case 'v'
        Release = vect_sum(Exit_2nm,20);
        %Release = vect_sum(Exit_10nm,20); 
        %Release = vect_sum(Exit_2nm_3,20); 
        Release = Release'; %Make Release a column vector
        initial =Release(1);
    otherwise
        OK = 0;
        disp('Enter V or I');
end
%%

%initial = 4000    $Set initial to a fixed amount if needed.
i = 60; j = 101; k = 101;  %Size of C and S matrices
%tin=input('Initial Timestep \n');
%tf=input('Final Timestep \n');
t = tin;

if tin >= tf
    disp('Initial timestep must be less than final');
    OK = 0;
elseif tin > 1 %add additional space to existing variables
    Ctemp = C; Stemp = S;
    %global S C
    C = Ctemp; S = Stemp;
    Cleft_Conc=cat(1,Cleft_Conc,zeros(tf-tin,1));
    Cleft_Zone=cat(1,Cleft_Zone,zeros(tf-tin,1));
    Cconc=cat(1,Cconc,zeros(tf-tin,1));
    R_dual=cat(1,R_dual,zeros(tf-tin,16));
else %start from scratch
    %global S C
    %pos = input('Release over which receptor, numeric value only \n');
    S = Smat(i,j,k,cleft_geo,cleft,zone);
    C = Cmat(i,j,k,initial,Release_loc);
    R_dual=zeros(tf,32);
    [i_n j_n k_n]=size(C);
    Ctemp=zeros(i_n,j_n,k_n);
    Cconc=zeros(tf+2,1);
    Cconc(2)=sum(sum(sum(C)));
    Cleft_Conc=zeros(tf+2,1);
    Cleft_Conc(2)=sum(sum(sum(C(29:31,21:81,21:81))));
    Cleft_Zone=zeros(tf+2,1);
    Cleft_Zone(2)=sum(sum(sum(C(29:31,34:67,34:67))));
end

%Cconc will be a function of t, that gives the total concentration of the
%3-dimensional array to check for conservation of glutamat.
Cconc(t+2)=sum(sum(sum(C)));
Cleft_Conc(t+2)=sum(sum(sum(C(29:31,21:81,21:81))));
Cleft_Zone(t+2)=sum(sum(sum(C(29:31,34:67,34:67))));

tic;


while OK == 1
    
%Start with the exterior
i=[2:i_n-1]; %same value for all for whole exterior region

%Region 1 and 2
k=[2:k_n-1]; j=[2:20,82:j_n-1]; 
    Ctemp(i,j,k) = diffALL(i,j,k,D1);
%Region 3 and 4
k=[2:20,82:k_n-1]; j=[21:81];
    Ctemp(i,j,k) = diffALL(i,j,k,D1);
    
%Interior, we first start with the 'outer' region of the cleft
%Boundary between exteriors and interior, the corners are dealt with
%separately.
i=29:31;    %i is constant for this region
j=21;  k=[22:80];
    Ctemp(i,j,k) = diffW(i,j,k,D2,Ddiff12); %West Boundary
j=81; 
    Ctemp(i,j,k) = diffE(i,j,k,D2,Ddiff12); %East Boundary
j=[22:80]; k=21;  
    Ctemp(i,j,k) = diffS(i,j,k,D2,Ddiff12); %South Boundary
k=81;
    Ctemp(i,j,k) = diffN(i,j,k,D2,Ddiff12); %North Boundary

%Corners with free space and cleft
j=21; k=21; 
    Ctemp(i,j,k) = diffSW(i,j,k,D2,Ddiff12); %SW Corner
j=81; k=21; 
    Ctemp(i,j,k) = diffSE(i,j,k,D2,Ddiff12); %SE Corner
j=21; k=81; 
    Ctemp(i,j,k) = diffNW(i,j,k,D2,Ddiff12); %NW Corner
j=81; k=81; 
    Ctemp(i,j,k) = diffNE(i,j,k,D2,Ddiff12); %NE Corner

%Outer area in cleft
%1st and 2nd outer region
j=[22:33,68:80]; k=[22:80];
    Ctemp(i,j,k) = diffALL(i,j,k,D2);
%3rd and 4th outer region    
j=[34:67];       k=[22:33,68:80];
    Ctemp(i,j,k) = diffALL(i,j,k,D2);

%Inner Area of cleft
j=[35:66]; k=[35:66];
    Ctemp(i,j,k) = diffALL(i,j,k,D3);
    
%Boundary between inner and outer area of synaptic cleft
j=34; k=[35:66];
    Ctemp(i,j,k) = diffW(i,j,k,D3,Ddiff23); %West boundary
j=67;  
    Ctemp(i,j,k) = diffE(i,j,k,D3,Ddiff23); %East boundary
j=35:66; k=34;  
    Ctemp(i,j,k) = diffS(i,j,k,D3,Ddiff23); %South boundary
k=67; 
    Ctemp(i,j,k) = diffN(i,j,k,D3,Ddiff23); %North boundary
    
%Now consider the four corners of the inner/outer cleft
j=34; k=34; 
    Ctemp(i,j,k) = diffSW(i,j,k,D3,Ddiff23); %SW Corner
j=67; k=34; 
    Ctemp(i,j,k) = diffSE(i,j,k,D3,Ddiff23); %SE Corner
j=34; k=67; 
    Ctemp(i,j,k) = diffNW(i,j,k,D3,Ddiff23); %NW Corner
j=67; k=67; 
    Ctemp(i,j,k) = diffNE(i,j,k,D3,Ddiff23); %NE Corner

%Updating C with Ctemp
C=Ctemp;
if synapse_size == 600
    R_dual(t+1,1)=C(31,26,76);R_dual(t+1,2)=C(31,38,76);R_dual(t+1,3)=C(31,51,76);R_dual(t+1,4)=C(31,64,76);R_dual(t+1,5)=C(31,76,76);
    R_dual(t+1,6)=C(31,26,64);R_dual(t+1,7)=C(31,38,64);R_dual(t+1,8)=C(31,51,64);R_dual(t+1,9)=C(31,64,64);R_dual(t+1,10)=C(31,76,64);
    R_dual(t+1,11)=C(31,26,51);R_dual(t+1,12)=C(31,38,51);R_dual(t+1,13)=C(31,51,51);R_dual(t+1,14)=C(31,64,51);R_dual(t+1,15)=C(31,76,51);
    R_dual(t+1,16)=C(31,26,38);R_dual(t+1,17)=C(31,38,38);R_dual(t+1,18)=C(31,51,38);R_dual(t+1,19)=C(31,64,38);R_dual(t+1,20)=C(31,76,38);
    R_dual(t+1,21)=C(31,26,26);R_dual(t+1,22)=C(31,38,26);R_dual(t+1,23)=C(31,51,26);R_dual(t+1,24)=C(31,64,26);R_dual(t+1,25)=C(31,76,26);
    
    
elseif synapse_size == 400
    R_dual(t+1,1)=C(31,25,57);R_dual(t+1,2)=C(31,33,57);R_dual(t+1,3)=C(31,41,57);R_dual(t+1,4)=C(31,49,57);R_dual(t+1,5)=C(31,57,57);
    R_dual(t+1,6)=C(31,25,49);R_dual(t+1,7)=C(31,33,49);R_dual(t+1,8)=C(31,41,49);R_dual(t+1,9)=C(31,49,49);R_dual(t+1,10)=C(31,57,49);
    R_dual(t+1,11)=C(31,25,41);R_dual(t+1,12)=C(31,33,41);R_dual(t+1,13)=C(31,41,41);R_dual(t+1,14)=C(31,49,41);R_dual(t+1,15)=C(31,57,41);
    R_dual(t+1,16)=C(31,25,33);R_dual(t+1,17)=C(31,33,33);R_dual(t+1,18)=C(31,41,33);R_dual(t+1,19)=C(31,49,33);R_dual(t+1,20)=C(31,57,33);
    R_dual(t+1,21)=C(31,25,25);R_dual(t+1,22)=C(31,33,25);R_dual(t+1,23)=C(31,41,25);R_dual(t+1,24)=C(31,49,25);R_dual(t+1,25)=C(31,57,25);
  
elseif synapse_size == 200
    R_dual(t+1,1)=C(31,23,39);R_dual(t+1,2)=C(31,27,39);R_dual(t+1,3)=C(31,31,39);R_dual(t+1,4)=C(31,35,39);R_dual(t+1,5)=C(31,39,39);
    R_dual(t+1,6)=C(31,23,35);R_dual(t+1,7)=C(31,27,35);R_dual(t+1,8)=C(31,31,35);R_dual(t+1,9)=C(31,35,35);R_dual(t+1,10)=C(31,39,35);
    R_dual(t+1,11)=C(31,23,31);R_dual(t+1,12)=C(31,27,31);R_dual(t+1,13)=C(31,31,31);R_dual(t+1,14)=C(31,35,31);R_dual(t+1,15)=C(31,39,31);
    R_dual(t+1,16)=C(31,23,27);R_dual(t+1,17)=C(31,27,27);R_dual(t+1,18)=C(31,31,27);R_dual(t+1,19)=C(31,35,27);R_dual(t+1,20)=C(31,39,27);
    R_dual(t+1,21)=C(31,23,23);R_dual(t+1,22)=C(31,27,23);R_dual(t+1,23)=C(31,31,23);R_dual(t+1,24)=C(31,35,23);R_dual(t+1,25)=C(31,39,23);
else    
end

        
%R_Loc(26,:)=[51,44];R_Loc(27,:)=[51,32];R_Loc(28,:)=[57,44];R_Loc(29,:)=[70,32];


%Cconc will be a function of t, that gives the total concentration of the
%3-dimensional array to check for conservation of glutamat.
Cconc(t+2)=sum(sum(sum(C)));
Cleft_Conc(t+2)=sum(sum(sum(C(29:31,21:81,21:81))));
Cleft_Zone(t+2)=sum(sum(sum(C(29:31,34:67,34:67))));

if mod(t,10000) == 0
    disp('Total Concentration is')
    disp(Cconc(t+2))
    disp('Current Cleft Concentration is')
    disp(Cleft_Conc(t+2))
    disp('At timestep (dt)')
    disp(t)
    disp('R_c or R_e Concentration is')
    disp(R_dual(t+1,5));
    etime=toc/60;
    disp('Elapsed time (minutes) is:')
    disp(etime)
    disp('.'),disp('.')
end

%Add release to cleft
t=t+1;
if initial ~= 4000
    if pos == 13 && t <= size(Release,1)
        if synapse_size == 600%%#ok<AND2>
            C(29,51,51) = C(29,51,51) + Release(t);
        elseif synapse_size == 400
            C(29,41,41) = C(29,41,41) + Release(t);
        elseif synapse_size == 200
            C(29,31,31) = C(29,31,31) + Release(t);
        end
    elseif pos == 25 && t <= size(Release,1)
        if synapse_size == 600%%#ok<AND2>
            C(29,76,26) = C(29,76,26) + Release(t);
        elseif synapse_size == 400
            C(29,57,25) = C(29,57,25) + Release(t);
        elseif synapse_size == 200
            C(29,39,23) = C(29,39,23) + Release(t);
        end
    end
end


if t == tf
    OK = 0;
end

end
conctime=toc/60;

% Total concentration plot
plot(Cconc)
set(gca,'XTickLabel',[0 0.02 0.04 0.06 0.08 0.1 0.12]);
% Concetrations on each receptor
figure
subplot(5,5,1);plot(R_dual(:,1));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,2);plot(R_dual(:,2));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,3);plot(R_dual(:,3));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,4);plot(R_dual(:,4));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,5);plot(R_dual(:,5));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,6);plot(R_dual(:,6));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,7);plot(R_dual(:,7));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,8);plot(R_dual(:,8));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,9);plot(R_dual(:,9));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,10);plot(R_dual(:,10));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,11);plot(R_dual(:,11));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,12);plot(R_dual(:,12));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,13);plot(R_dual(:,13));title('R_{center}');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,14);plot(R_dual(:,14));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,15);plot(R_dual(:,15));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,16);plot(R_dual(:,16));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,17);plot(R_dual(:,17));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,18);plot(R_dual(:,18));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,19);plot(R_dual(:,19));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,20);plot(R_dual(:,20));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,21);plot(R_dual(:,21));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,22);plot(R_dual(:,22));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,23);plot(R_dual(:,23));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,24);plot(R_dual(:,24));title('');set(gca,'XTickLabel',[0 0.05 0.1]);
subplot(5,5,25);plot(R_dual(:,25));title('R_{edge}');set(gca,'XTickLabel',[0 0.05 0.1]);






