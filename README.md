Modeling-of-Neurotransmissions for Indepentdent Evoked and Spontaneous Release to Receptors

1. First step : Concentrations of neurotransmitters after release from presynaptic neurons(Diffusion process)

  Setting up the parameters for diffusion of neutrotransmitters(e.g. Glutamate molecues)
  
  Run 
  
  

Frist 

%% Opening probability procedure;

%% for Center release(R6)
[Time_R6,probM_R6,probL_R6,totM_R6,totL_R6]=NMDA_kinetics(R_dual(:,6),Ndt,av,inc);
[Time_R10,probM_R10,probL_R10,totM_R10,totL_R10]=NMDA_kinetics(R_dual(:,10),Ndt,av,inc);
[Time_R14,probM_R14,probL_R14,totM_R14,totL_R14]=NMDA_kinetics(R_dual(:,14),Ndt,av,inc);

% Center release
figure
plot(probM_6,'r-');hold on;
plot(probM_16,'b-.');hold on;
figure(gcf)
set(gca,'XTickLabel',[0 10 20 30 40 50 60 70 80 90 100]);
title('Open probabilities released on center site')
xlabel('Time')
ylabel('Open Probability')

% For measurement
PeakOP_center_600_1=zeros(3,1);
PeakOP_center_600_1(1,1) = max(probM_R6)/max(probM_R6);
PeakOP_center_600_1(2,1) = max(probM_R10)/max(probM_R6);
PeakOP_center_600_1(3,1) = max(probM_R14)/max(probM_R6);


% For measurement
PeakOP_center_600_2=zeros(3,1);
PeakOP_center_600_2(1,1) = max(probM_R6)/max(probM_R6);
PeakOP_center_600_2(2,1) = max(probM_R11)/max(probM_R6);
PeakOP_center_600_2(3,1) = max(probM_R16)/max(probM_R6);


figure
%plot(X_min_c_600,Y_min_c_600,'or',t_min_c_600,pchip(X_min_c_600,Y_min_c_600,t_min_c_600),'r-');hold on;
plot(PeakOP_center_600','or',pchip(X_max_c_600,PeakOP_center_600',t_max_c_600),'r-','MarkerSize',10);hold on;
%plot(X_min_c_400,Y_min_c_400,'or',t_min_c_400,pchip(X_min_c_400,Y_min_c_400,t_min_c_400),'b--');hold on;
plot(X_max_c_600,PeakOP_center_600_highaffiC','*b',t_max_c_600,pchip(X_max_c_600,PeakOP_center_600_highaffiC',t_max_c_600),'b--','MarkerSize',10);hold on;
%plot(X_min_c_200,Y_min_c_200,'or',t_min_c_200,pchip(X_min_c_200,Y_min_c_200,t_min_c_200),'m-.');hold on;
plot(X_max_c_600,PeakOP_center_600_highaffiE','^m',t_max_c_600,pchip(X_max_c_600,PeakOP_center_600_highaffiE',t_max_c_600),'m-.','MarkerSize',10)
ylim([0,1.05]);
title('Ratio of Peak Open Probabilities for Center Release with Diffusion Constraints')
ylabel('Peak OP at each location/Peak OP over center site')
xlabel('Distance(nm)')
legend('Base Model','(600nm)','High Affinity','Center','High Affinity','Edge')
grid;





%% Edge release (R16)
[Time_R16,probM_R16,probL_R16,totM_R16,totL_R16]=NMDA_kinetics(R_dual(:,16),Ndt,av,inc);
[Time_R11,probM_R11,probL_R11,totM_R11,totL_R11]=NMDA_kinetics(R_dual(:,11),Ndt,av,inc);
[Time_R6,probM_R6,probL_R6,totM_R6,totL_R6]=NMDA_kinetics(R_dual(:,6),Ndt,av,inc);
[Time_R1,probM_R1,probL_R1,totM_R1,totL_R1]=NMDA_kinetics(R_dual(:,1),Ndt,av,inc);


%Edge release
figure
plot(probM_6,'r-');hold on;
plot(probM_16,'b-.');hold on;
figure(gcf)
set(gca,'XTickLabel',[0 10 20 30 40 50 60 70 80 90 100]);
title('Open probabilities released on edge site')
xlabel('Time')
ylabel('Open Probability')


% For measurement
PeakOP_edge_600=zeros(4,1);
PeakOP_edge_600(1,1) = max(probM_R16)/max(probM_R16);
PeakOP_edge_600(2,1) = max(probM_R11)/max(probM_R16);
PeakOP_edge_600(3,1) = max(probM_R6)/max(probM_R16);
PeakOP_edge_600(4,1) = max(probM_R1)/max(probM_R16);


pchip(PeakOP_center_600)


      x = 0:3;
      y = [0 0.2 0.4 0.6 0.8 1];
      t = 0:.01:3;
      
      plot(PeakOP_center_600,'or',pchip(x,PeakOP_center_600,t),'r-','MarkerSize',10);hold on;

      plot(x,y,'o',t,pchip(PeakOP_center_600), pchip(PeakOP_center_600)])
      legend('data','pchip','spline', 'Location', 'Best')
