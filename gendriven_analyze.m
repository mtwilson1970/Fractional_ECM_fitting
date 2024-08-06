%
% 
%
%gendriven_analyze.m
%
% Marcus Wilson. Fit a CPE_CPE_R model to V(t) and I(t) data, by different
% days. 
%
%24 June 2024.  Allow for charge and discharge Resistances to be different
%if necessary.
%p9373

close all;
clear;

format long;   tic;

%%%%%%%%%Set up parameters%%%%%%%%%%%%%%%%%

%%%%More details are defined in function fit_CPE_CPE_R_as_func, including
%%%%the data file we want to use. %%%%%%%%%%%%%%%%

backdays=0;   %How many days of history do we take into account with fitting? Default = 0 which means each day is treated independently. 
day_range=[1 2]; %range of days to fit [first last]. 2 days will do the Grydrive4.mat example.
%96 days for U1drive3.mat
%Choose to use (say) 10-day periods rather than one day (as defined in
%fit_CPE_CPE_R_as_func.m) we will divide day range by ten - i.e. [1 9] will
%do nine lots of ten days. 



alphas1=[0.97:0.005:1.00];   %sensible ones. Range of alpha values considered, [min:stepsize:max]
alphas2=[0.12:0.01:0.30];  

downing=20;   %How much do we downsample data by to keep memory requirements reasonable?

%%%%%%%%%%%%%%%%%%%%%End set up parameters%%%%%%%%%%%%%%%%%%%%%



%Vance's impedance data from electrical impedance spectroscopy

%for U1drive
ffEIS = [1e-5; 2e-5; 5e-5; 1e-4; 2e-4; 5e-4; 1e-3; 2e-3; 5e-3; 1e-2; 2e-2; 5e-2; 1e-1];
%modZEIS = [1.1; 0.64; 0.28; 0.17; 0.10; 0.07; 0.06; 0.053; 0.047; 0.043; 0.041; 0.039; 0.038];
%phaseEIS = -[83; 79; 73; 63; 48; 30; 22; 17; 14; 11; 8; 5; 4];

%Actual values taken from machine not graph
modZEIS=[1.195; 0.644; 0.283; 0.160; 0.101; 0.0705; 0.0608; 0.0555; 0.0480; 0.0437; 0.0411; 0.0391; 0.0381];
phaseEIS=-[81.7; 79.4; 73.2; 62.40; 48.0; 30.1; 21.5; 16.5; 13.4; 10.3; 7.44; 4.89; 3.70];
modZEIS2= [1.213; 0.652; 0.290; 0.167; 0.110; 0.0794; 0.0688; 0.0622; 0.0537; 0.0488; 0.0460; 0.0438; 0.0427];
phaseEIS2= -[81.1; 78.1; 71.0; 59.6; 45.2; 28.8; 20.9; 16.6; 13.5; 10.3; 7.42; 4.89; 3.68];


ff=logspace(-6,1,22);

%Set up blank arrays that will be filled in. 
Rch_fit = NaN*zeros(1,max(day_range(2),1));  %charge
    Rdi_fit = Rch_fit; %discharge
    R_fit = Rch_fit;    %mean
    a1_fit = Rch_fit;
    C1_fit = Rch_fit;
    a2_fit = Rch_fit;
    C2_fit = Rch_fit;
    V0_fit = Rch_fit;
    lsq_fit=Rch_fit;
  Zmat = NaN*zeros(11,length(ff));
  
  
  
  
 for day=day_range(1):day_range(2);    %Loop over all days selected
     
     disp(['Evaluating for day ' num2str(day)])
    
    %go over all days
    [lsq,Rfittedch,Rfitteddi,alpha_best_1,CFfitted1,alpha_best_2,CFfitted2,constV]=fit_CPE_CPE_R_as_func(day,backdays,downing,alphas1,alphas2,true)
    
    if (day==0)
        %if we are fitting all days, increment it so we have an index
        %bigger than one
        dayp=1;
    else
        dayp=day;  
    end
    
    %load up best fits
    Rch_fit(dayp) = Rfittedch;
    Rdi_fit(dayp) = Rfitteddi;
    a1_fit(dayp) = alpha_best_1;
    C1_fit(dayp) = CFfitted1;
    a2_fit(dayp) = alpha_best_2;
    C2_fit(dayp) = CFfitted2;
    V0_fit(dayp) = constV;
    lsq_fit(dayp) = lsq;
    
    R_fit(dayp) = 0.5*(Rfittedch + Rfitteddi);   %mean R
    
    %find impedance
    %Not strictly correct for different ch and di currents. Then impedance
    %doesn't really have a meaning
    Z = 0.5*(Rfittedch+Rfitteddi) + (CFfitted1 * (2*pi*ff*j).^alpha_best_1).^(-1) ...
        + (CFfitted2 * (2*pi*ff*j).^alpha_best_2).^(-1);
    
    Z_mat(dayp,:) = Z(:);   %load up matrix.
    
    
 

    
    
    
 end

 %analyze
 
 save data_dump_long
 
 
 modZ = real( (Z_mat.*conj(Z_mat)) ).^(0.5); %Modulus of Z
 phaseZ = (180/pi)*(atan2(imag(Z_mat),real(Z_mat)));  %put it in degrees. Phase of Z
 
 figure(9);   subplot(2,1,1)
 
 loglog(ffEIS,modZEIS,'bo'); grid on; hold on;  %plot Vance's EIS result for U1
 
 %16 jan 2024. Do the next Z one from 12 January 2024.
  loglog(ffEIS,modZEIS2,'kx'); grid on; hold on;  
 
  
 xlabel('frequency (Hz)')
 ylabel('|Z| (ohm)')
 subplot(2,1,2)
 semilogx(ffEIS,phaseEIS,'bo'); grid on; hold on;  %Vance's phase
 semilogx(ffEIS,phaseEIS2,'kx'); grid on; hold on;  %2nd lot of EIS on U1 drive
 
 xlabel('frequency (Hz)')
 ylabel('phase of Z (degrees)')
 
 for day=1:max(day_range(2),1);
     cv=day/day_range(2);
    subplot(2,1,1); loglog(ff,modZ(day,:),'-','color',[1-cv 4*cv*(1-cv) cv]); grid on; hold on;
    subplot(2,1,2); semilogx(ff,phaseZ(day,:),'-','color',[1-cv 4*cv*(1-cv) cv]); grid on; hold on;
     
     %find Z at knee freq
     knee_omega = (R_fit(day)*C1_fit(day))^(-1/a1_fit(day));  
     Z_at_knee(day) =  R_fit(day) + (C1_fit(day) * (knee_omega*j)^a1_fit(day))^(-1) ...
        + (C2_fit(day) * (knee_omega*j)^a2_fit(day))^(-1);

     
 end
 if (day_range(2) == 0)
     %Just the fit to the whole sequqnce
     legend('from EIS', 'all days');
 else
 legend('from EIS','day 1','day 2','day 3','day 4','day 5','day 6',...
     'day 7', 'day 8', 'day 9', 'day 10', 'day 11');
 end
 %Record the best paras   
 subplot(2,1,1);
 title(['alpha1=' num2str(mean(a1_fit))  ...
        ', CF1=' num2str(mean(C1_fit)) ' A s^{a}/V' ...
        ', alpha2=' num2str(mean(a2_fit)) ...
        ', CF2=' num2str(mean(C2_fit)) ' A s^{a}/V' ...
        ', Rch=' num2str(mean(Rch_fit)) ' ohm' ...
        ', Rdi=' num2str(mean(Rdi_fit)) ' ohm' ...
        ', Vc=' num2str(mean(V0_fit)) ' V'])
 
 
%  %Relative fits
%  ave_by_log = exp(mean(log(modZ),1));
%  ave_by_log_matrix = ones(11,1)*ave_by_log;
%  rel_Z = modZ./ave_by_log_matrix;
%  
%  figure(10)
%  
%  for day=1:max(day_range(2),1)
%      cv=day/day_range(2);
%      semilogx(ff,rel_Z(day,:),'-','color',[ 1-cv 4*cv*(1-cv) cv]); grid on; hold on;
%  end
%  legend('day 1','day 2','day 3','day 4','day 5','day 6',...
%      'day 7', 'day 8', 'day 9', 'day 10', 'day 11');
%  xlabel('frequency (Hz)')
%  ylabel('|Z| relative to mean (by log) across days')
 
figure(11)
plot([1:max(day_range(2),1)],sqrt(real(Z_at_knee.*conj(Z_at_knee))),'k-')
xlabel('day')
ylabel('mod impedance at knee (ohm)')

%Find cross over points.  8 December 2022



f1=(1/(2*pi))*(C1_fit./C2_fit).^( (a2_fit - a1_fit).^(-1) );
f2=(1/(2*pi))*(C2_fit.*R_fit).^(-(a2_fit.^(-1))); 

figure(12); grid on; hold on
plot([1:max(day_range(2),1)],log10(f1),'kx-'); 
plot([1:max(day_range(2),1)],log10(f2),'bo-');
xlabel('day')
ylabel('log10(crossover frequencies (Hz))')

%vance's total gendrive one
%Rv=0.027; a1v=0.986; CF1v=7573; a2v=0.055; CF2v=8.0;

%vance's EIS one
%Rv=0.161; a1v=0.960; CF1v = 5419; a2v=0.190; CF2v=94;

%vance's cycle one
%Rv=0; a1v=0.962; CF1v=6223; a2v=0.018; CF2v=5.8;

%my ave  (average values across fitting all days)
Rv=0.1521; a1v=0.9913; CF1v=7986; a2v=0.2082; CF2v=72.6;

f1v=(1/(2*pi))*(CF1v/CF2v)^(1/(a2v-a1v));
f2v=(1/(2*pi))*(CF2v*Rv)^(-1/a2v); 
plot([1:max(day_range(2),1)],ones(1,max(day_range(2),1))*log10(f1v),'k--');
plot([1:max(day_range(2),1)],ones(1,max(day_range(2),1))*log10(f2v),'b--');
legend('f1','f2','f1(ave day)','f2(ave day)')

toc
