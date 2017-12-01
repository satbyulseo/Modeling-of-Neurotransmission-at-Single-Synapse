Modeling-of-Neurotransmissions for Indepentdent Evoked and Spontaneous Release to Receptors

1. Concentrations of neurotransmitters after release from presynaptic neurons(Diffusion process)

  Run Difussion_Release.m
  
    1)  Size of synapse, 200, 400, or 600 nm? %This is for managing size of synapse.
         Input:  numeric value either 600, 400, or 200 
    
    2)  Release over which receptor, numeric value only. 
         Input : 13 or 25
         %13 represents release of neurotransmitters from the center preynaptic terminal.
         %25 represnets release of neurotransmitters from the edge presynaptic terminal.
    
    3) Dglut value for outer region of cleft? 
        input : numeric value from 0.1-0.4 
        
       %This Diffusion constant represents the condition(Speed, barrier, viscosity, etc.) inside of cleft. 
    
    4) Dglut value for inner region of cleft? 
       Same as #3)
       
    5) Instantaneous, I, or Vesicle release, V?
       Input: either I or V 
       
       % I represents normal release process
       % V represents release process with limited amount of neurotransmitter release. 
       
 
2. Opening probability for NMDA receptor using Kinetic model.

    1) Setting for parameter "C" using output from diffusion process
       
       C = R_dual(:,'Numeric value')
       
       %"Numeric Value" is raning from 1 to 25 represeting the location on the postsynaptic.
       %You can input any numbers that you want to look at the opening probability.
       
       Run 
       
      ex) [Time_13,probM_13,totM_13]=NMDA_kinetics(R_dual(:,13),Ndt,av,inc);
          [Time_25,probM_25,totM_25]=NMDA_kinetics(R_dual(:,25),Ndt,av,inc);
       

    
    3) Plot
    
     % Center release
        
        figure
        plot(probM_13,'r-');hold on;
        plot(probM_25,'b-.');hold on;
        figure(gcf)
        set(gca,'XTickLabel',[0 10 20 30 40 50 60 70 80 90 100]);
        title('Open probabilities released From Center Presynaptic Terminal');
        xlabel('Time');
        ylabel('Open Probability');



