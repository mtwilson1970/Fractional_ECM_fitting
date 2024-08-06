function [lsq,Rfittedch,Rfitteddi,alpha_best_1,CFfitted1,alpha_best_2,CFfitted2,constV]=fit_CPE_CPE_R_as_func(whichday,backdays,downing,alphas1,alphas2,do_we_plot)

%As a function, gets the best fit paras
%
% Marcus Wilson
% 10 November 2022
% Fitting a CPE-CPE-R model to Vance's gendrive data.
% Development on previous fit_CPE_R code.
% 
% whichday gives us which of the days (1-11). 0 is all days
% backdays gives how many previous days data do we take into account?
% downing is downsampling rate (about 50)
% alphas1 gives range of possible alpha 1.
% alphas2 gives range of possible alpha 2
% do_we_plot is a flag - true and we plot some subsidiary graphs (and it
% takes more time), false and we don't.

%24 June 2024. Have Rfittedch and Rfitteddi options, so can have two
%different R values if we choose so. 

close all; 


%%%Input parameters and options defined here%%%%%%%%%%%%%%%
one_or_two=2;   %Do we fit one or two CPEs? %by default do 2 CPEs unless overridden (e.g. U1drive)
%A value of 3 fits 2 CPEs and also a different Rch and Rdi value.   Use 2
%unless for other research. 

%Which data sequence do we analyze?   U1drive is done for Wilson et al
%(2024) manuscript "Early detection of Li-ion cell failure from
%current-voltage time-series data. 

use_GenDriveN=false;
use_tc=false;
use_cyc=false;
use_U1drive=false;    %Use this one for 2024 manuscript Wilson et al. 
use_Grydrive=true;    %Use this one as an example

%%Now look for variable_name in the relevant section below, for the name of
%%the load file (e.g. U1drive3)

%Also look for day_stride which tells us how many days get grouped
%together. (Default=1, for one day at a time). 

fit_constant_voltage=true;   %do we fit the constant voltage (true) or keep it as initial value (false)

%%%%%%%%%%%%End set up parameters%%%%%%%%%%%%%%%%%%%

%Could use these
%file_to_use_list={'tc133xxx_15mar22' 'tc153xxx_15mar22' 'tc173xxx_15mar22' 'tc193xxx_15mar22'}


if (use_GenDriveN)

load GenDriveN.mat;
%load Mauve4SHCRC3.mat

if (whichday == 0)
    alldays=true;
else
alldays=false;   %do we do all days. If not, whichday is important
end



if (alldays)
    x=99;
    y=1;
else
     x=100/(backdays+1)-1;  %we assess the fit to the last x% of the data sequence
    y=1;  %but not the last y%  
    
  
  %  x=15.5; y=1;
    
end






%resample onto fixed time points
%16 Jan 2023. Keep the original points in place. 

deltat_from_file = min(tt(2:end)-tt(1:end-1));


Fs = 1/deltat_from_file;       %sample at this rate.
[Is,tts]=resample(I,tt,Fs);
[Vs,tts]=resample(V,tt,Fs);   





%10 November 2022. Now use downsample not resample
%Ir=resample(Is,1,downing);
%Vr=resample(Vs,1,downing);
%ttr=resample(tts,1,downing);

Ir=downsample(Is,downing);
Vr=downsample(Vs,downing);
ttr=downsample(tts,downing);

%Now we can clear the big variables

clear GenDrive I V tt; 
%Deltat
deltat=ttr(2)-ttr(1);   %need to assume it's constant.   Should be if resampled.

%index for start of a use. Battery is 'rested' at zero current before these
%indices





if (alldays)
    lengthmin=1;
    lengthmax=length(Vr);
else
    index_starts=(50/downing)*[20e3 40e3 60e3 80e3 100e3 120e3 140e3 161e3 183e3 204e3 224e3 241e3]-1;

    lengthmin=index_starts(whichday-backdays);   %start val   Which indices do we take
    lengthmax=index_starts(whichday+1)-1;    %end val
end

% %no chopping up into days
% lengthmin=1
% lengthmax=length(Vr);

Ir=Ir(lengthmin:lengthmax); Vr=Vr(lengthmin:lengthmax); ttr=ttr(lengthmin:lengthmax);


end


if (use_cyc)
    %Use the Cyc,Shal,Mysa etc series of data to look at the Cycs.
    
    disp('***Using cyc data***')
    
    if (whichday==1)
        load Cyc310xx_23apr22.mat
        basename='Cyc310xx';
    end
    if (whichday==2)
        load Cyc400xx_28apr22.mat
        basename='Cyc400xx';
    end
    if (whichday==3)
        load Cyc430xx_3may22.mat
        basename='Cyc430xx';
    end
    if (whichday==4)
        load Cyc520xx_10may22.mat
        basename='Cyc520xx';
    end
    if (whichday==5)
        load Cyc580xx_17may22.mat
        basename='Cyc580xx';
    end
     if (whichday==6)
        load Cyc760xx_28may22.mat
        basename='Cyc760xx';
    end
    
    %resample the data here 
    tt=eval([basename '(:,1)']);
    V=eval([basename '(:,2)']);
    I=eval([basename '(:,3)']);
    
    [Is,tts]=resample(I,tt);   %make equal deltat
    [Vs,tts]=resample(V,tt);  

    Ir=downsample(Is,downing);
    Vr=downsample(Vs,downing);
    ttr=downsample(tts,downing);
    
    deltat=ttr(2)-ttr(1);
    x=99; y=1;
end


if (use_tc)
  load tc133xxx_15mar22.mat
  
  %resample onto fixed time points
  tt=tc133xxx(:,1);
  V=tc133xxx(:,2);
  I=tc133xxx(:,3);
  
[Is,tts]=resample(I,tt);   %make equal deltat
[Vs,tts]=resample(V,tt);  
downing=100;   %resample
Ir=downsample(Is,downing);
Vr=downsample(Vs,downing);
ttr=downsample(tts,downing);

x=99;  %we fit to the last x% of the data sequence
y=1;  %but not the last y%


  clear tc133xxx_15mar22.mat; 
   %Deltat
   deltat=ttr(2)-ttr(1);   %need to assume it's constant.   Should be if resampled.
end
  

if (use_U1drive)
    %7 November 2023
    %Using the U1drive data sequence
    
    day_stride=1;   %do in 10 day incs not 1 day incs. Then day_stride=10;
    
   % loadname='U1comb';  %(what the .mat file is called, without the .mat)
   % variablename='U1comb'  %(what the variable is called

    loadname='U1drive3';   %Data in this case is the U1drive3.mat file
    variablename='U1drive3';    %Variable is called U1drive3
    
    load(loadname);
    eval(['tt=' variablename '(:,1);']);
    eval(['V=' variablename '(:,2);']);
    eval(['I=' variablename '(:,3);']);

    if (length(tt) > 69940000)
        tt(69940543)=(tt(69940542) + tt(69940544))/2;  %p9304 correct dodgy point
        tt(69940061)=(tt(69940060) + tt(69940062))/2;  % a nan
    end

    deltat=mean(tt(5:end-5)-tt(4:end-6))

    day_length=86400/deltat  %how many samples in a day

    if (whichday(1)==0)
        %All sequence
        dayst=1; dayen=length(V)-1;   %ignore last point, sometimes a NaN;
    else
        dayst=round(  (whichday-1)*day_stride*day_length+1 ) %starting point for day
      %  dayen=day(2)*day_length; %end point for day
         dayen=round(dayst+day_stride*day_length-1)  %index for end point of single day
    end

    Ir=downsample(I(dayst:dayen),downing); 
    Vr=downsample(V(dayst:dayen),downing);  %construct I and V sequences
    ttr=downsample(tt(dayst:dayen),downing);

    samplef=1/(deltat*downing)   %frequency samples
     deltat=ttr(2)-ttr(1)   %need to assume it's constant.   Should be if resampled. Needs to be the resampled one
    x=99;   %fit to last x% but not to last y%
    y=1;

    clear U1comb
     
    one_or_two=2;    %3 for charge and discharge resistance separately
end

if (use_Grydrive)
    %7 November 2023
    %Using the Grydrive data sequence
    
    day_stride=1;   %do in 10 day incs not 1 day incs. Then day_stride=10;
    
    loadname='Grydrive4';  %(what the .mat file is called, without the .mat)
    variablename='Grydrive4'  %(what the variable is called


    load(loadname);
    eval(['tt=' variablename '(:,1);']);
    eval(['V=' variablename '(:,2);']);
    eval(['I=' variablename '(:,3);']);

   % %Multiply current by 1000 to put it in milliamps. Helps matrix calcs?
   % I=I*1000;

    deltat=mean(tt(5:end-5)-tt(4:end-6))

    day_length=86400/deltat  %how many samples in a day

    if (whichday(1)==0)
        %All sequence
        dayst=1; dayen=length(V)-1;   %ignore last point, sometimes a NaN;
    else
        dayst=round(  (whichday-1)*day_stride*day_length+1 ); %starting point for day
      %  dayen=day(2)*day_length; %end point for day
         dayen=round(dayst+day_stride*day_length-1);  %index for end point of single day
    end

    %Correct.  Delay by 0.6592 days
    %delay_by_days=0.6592;    %How many days to pause by (see p9778)
    delay_by_days=0;
    delay_by_index = delay_by_days*86400/deltat;    %How many indices to pause by?
    dayst=dayst+delay_by_index;
    dayen=dayen+delay_by_index;  
    
    
    Ir=downsample(I(dayst:dayen),downing); 
    Vr=downsample(V(dayst:dayen),downing);  %construct I and V sequences
    ttr=downsample(tt(dayst:dayen),downing);

    samplef=1/(deltat*downing)   %frequency samples
     deltat=ttr(2)-ttr(1)   %need to assume it's constant.   Should be if resampled. Needs to be the resampled one
    x=99;   
    y=1;
    clear Grydrive
     
    one_or_two=2;
end


%For each alpha, we'll find a best fit R and CF. Then which alpha minimizes

CFstandard=1; Rstandard=1;   %standard paras to call

best_index=0;
smallest_so_far=99999999;   %a big number



lsq_array=NaN*zeros(length(alphas1),length(alphas2)); %initialize an array of results


get_V_matrix_1=NaN*zeros( length(Vr),length(alphas1) );  %initialize a matrix
get_V_matrix_2=NaN*zeros( length(Vr),length(alphas2) );  %initialize a matrix


for j=1:length(alphas1)
    
    alpha=alphas1(j) ;
    
    disp(['Evaluating for alpha=' num2str(alpha)])
    
    
    %Find the fraccap response
    
   %Find the fraccap response
    getV=find_V_with_alpha_CF(alpha,CFstandard,deltat,Ir);
    
    %Load up a matrix with this result in.
    get_V_matrix_1(:,j) = getV;

end
    %Now we have a matrix with the getV function for each alpha
    
    %Do second CPE
  for j=1:length(alphas2)
    
    alpha=alphas2(j); 
    
    disp(['Evaluating for alpha=' num2str(alpha)])
    
    
    %Find the fraccap response
    getV=find_V_with_alpha_CF(alpha,CFstandard,deltat,Ir);
    
    %Load up a matrix with this result in.
    get_V_matrix_2(:,j) = getV;

end  
    
    %Now find the Rs and CF that best fits with a pair of alphas
 suppress=~do_we_plot;    %don't display output - too  much of it
 
  alpha_best_1 = alphas1(1);
  alpha_best_2 = alphas2(1);   %select first ones as the best so far. 
  best_index_1 =1; best_index_2 = 1;
 
 for i=1:length(alphas1);
     
     for j=1:length(alphas2);
         
         
         
         disp(' ')
         disp(' ')
         alpha1=alphas1(i);
         alpha2=alphas2(j);
         disp(['Fitting for alpha1=' num2str(alpha1) ' and alpha2=' num2str(alpha2)]);
         
         
         getV1=get_V_matrix_1(:,i);   %get the first function
         getV2=get_V_matrix_2(:,j);   %get the second function
 
 
 [lsq,Rfittedch,Rfitteddi,CFfitted1,CFfitted2,constV,Vout]=CF_and_R_fit_twice(x,y,getV1,getV2,Ir,Vr,suppress,fit_constant_voltage,one_or_two);
    
 %8 Novemberr 2023. Reject anything if non physical - i.e. negative values,
 %by setting lsq to something big.
 if (Rfittedch < 0)
 %    lsq=99999999;
 end
 if (CFfitted1 < 0)
  %   lsq=99999999;
 end
 if (CFfitted2 < 0)
  %   lsq=99999999;
 end 
 if (constV < 0)
  %   lsq=99999999;
 end
 
 
    lsq_matrix(i,j) = lsq;
    if (lsq<smallest_so_far)
        %update new
        smallest_so_far = lsq;
        best_index_1=i;
        best_index_2=j;
        alpha_best_1=alphas1(i); Rfittedch_best=Rfittedch;    Rfitteddi_best = Rfitteddi;  CF_best_1 = CFfitted1;
        alpha_best_2=alphas2(j);  CF_best_2 = CFfitted2;

    end

    lsq_array(i,j)=lsq;
    
     end    %end second alpha loop
 end    %end first alpha loop

disp(' ')
disp(' ')
disp(['Best fitted alpha1 and alpha2 were ' num2str(alpha_best_1) ' and ' ...
    num2str(alpha_best_2)]);



%Now recalc the best case
suppress=false;   %don't suppress the output
getV1=get_V_matrix_1(:,best_index_1);   %get the first function
getV2=get_V_matrix_2(:,best_index_2);   %get the second function
 
[lsq,Rfittedch,Rfitteddi,CFfitted1,CFfitted2,constV,Vout]=CF_and_R_fit_twice(x,y,getV1,getV2,Ir,Vr,suppress,fit_constant_voltage,one_or_two);

if (do_we_plot)

figure
plot(ttr,Ir,'k-'); grid on;
xlabel('time (s)')
ylabel('actual current (A)')

figure
plot(ttr,Vr,'k-'); grid on; hold on;  %actual voltage
xlabel('time (s)')
ylabel('voltage (V)')
plot(ttr,Vout,'r-')  %modelled voltage
legend('actual', 'fitted')




disp(['rms fit is ' num2str(1000*sqrt(lsq)) ' mV']);

figure
contourf(alphas1,alphas2,1000*lsq_array'.^0.5,'kx-')
xlabel('alpha1 value')
ylabel('alpha2 value')
hold on;
plot3(alpha_best_1,alpha_best_2,100,'yx','MarkerSize',10)
title('r.m.s. fitting error (mV)')
colorbar

figure
plot(ttr,Vout-Vr,'b-')
xlabel('time (s)')
ylabel('modelled minus actual (V)')
grid on; hold on

end


end