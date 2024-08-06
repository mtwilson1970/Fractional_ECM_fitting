function [lsq,Rfittedch,Rfitteddi,CFfitted1,CFfitted2,constV,Vout]=CF_and_R_fit_twice(x,y,getV1,getV2,actualI,actualV,suppress,fit_const,one_or_two_CPE)
%
% Marcus Wilson
%
%10 November 2022   Fits two CPEs and an R.  Based on CF_and_R_fit


% This function fits a value of series resistor R and frac capacitor value
% CF given a standard VR for R=1  (i.e. actualI) and a standard VC for CF=1 (the
% getV.
%
% Inputs are getV (the standard voltage for a frac cap for a given alpha
% and for CF=1 unit, actualI (the actual current, A), and actualV (the
% actual voltage (V).
% suppress is true if we want to hide the text output
% fit_const is true if we fit the constant voltage part. If it's false we
% use the initial value for the voltage for this parameter.
% one_or_two_CPE = 1 if we fit just 1 CPE or 2 if we fit two. If it's not
% two we will fit 1.
%
% Outputs are lsq (the average lsq value, ave over number of points),
% Rfitted (the fitted resistor, ohms) and
% CF1 (the fitted frac cap, A s^alpha1 V-1) for first CPE of alpha1
% CF2 (the fitted frac cap) for second CPE of alpha2
% and Vout, the fitted voltage trace

%p8049. Use regression to fit V = a VC + b VR to a actualV

%fit to last x% of data but not to last y%

%p9373 24 June 2024. Allow fitting of two different R values for charge and
%discharge






N=length(getV1);
startpoint=ceil(N-(x/100)*N);
endpoint=ceil(N-(y/100)*N);

exclude_by_voltage=false;

if (exclude_by_voltage)
%%14 December 2023. Weight function
%%Don't count if actualV > 4.2 V
ww=((actualV<4.1) & (actualV>3.7)); 
%sumww=sum(ww);   disp(['sum of accepted terms is ' num2str(sumww)])

else

%exclude by charge
    
%Logan Cowie's idea 14 Dec 23 is to take out the top few points on coulomb
%count rather than voltage
Q=cumtrapz(actualI);
Qtop=max(Q); Qbottom=min(Q);   %max and min
Qn = (Q - Qbottom)/(Qtop - Qbottom);  %normalized between 0 (low) and 1 (max)
%ww=((Qn < 0.75) & (Qn > 0.25));   %take out topmost ones / bottom most ones
ww=(Qn<1.1);  %all of them


sumww=sum(ww);

end

%create terms


%Notation is V1 for first alpha sequence, V2 for second, V3 for I sequence
%24 June 2024.   Two I sequences now.   V3 is Icharge, V4 is Idischarge

%Ich and Idi
Ich = max(actualI,0);   %the positive ones
Idi = min(actualI,0);   %the negative ones



sumV1sq = sum(ww(startpoint:endpoint).*getV1(startpoint:endpoint).^2);
sumV2sq=  sum(ww(startpoint:endpoint).*getV2(startpoint:endpoint).^2);    %for second alpha
sumV3sq = sum(ww(startpoint:endpoint).*actualI(startpoint:endpoint).^2);  %only one R
sumV3chsq = sum(ww(startpoint:endpoint).*Ich(startpoint:endpoint).^2);  %charge
sumV4disq = sum(ww(startpoint:endpoint).*Idi(startpoint:endpoint).^2);  %discharge


sumV1V2 = sum(ww(startpoint:endpoint).*getV1(startpoint:endpoint).*getV2(startpoint:endpoint));
sumV1V3 = sum(ww(startpoint:endpoint).*getV1(startpoint:endpoint).*actualI(startpoint:endpoint));   %for only one R
sumV1V3ch = sum(ww(startpoint:endpoint).*getV1(startpoint:endpoint).*Ich(startpoint:endpoint));    %for charge
sumV1V4di = sum(ww(startpoint:endpoint).*getV1(startpoint:endpoint).*Idi(startpoint:endpoint));   %for discharge


sumV2V3 = sum(ww(startpoint:endpoint).*getV2(startpoint:endpoint).*actualI(startpoint:endpoint));  %for only one R
sumV2V3ch = sum(ww(startpoint:endpoint).*getV2(startpoint:endpoint).*Ich(startpoint:endpoint));   %for charge
sumV2V4di = sum(ww(startpoint:endpoint).*getV2(startpoint:endpoint).*Idi(startpoint:endpoint));  %for discharge

sumV3chV4di = sum(ww(startpoint:endpoint).*Ich(startpoint:endpoint).*Idi(startpoint:endpoint));   %should be zero because one or other will be



sumV1 = sum(ww(startpoint:endpoint).*getV1(startpoint:endpoint));
sumV2 = sum(ww(startpoint:endpoint).*getV2(startpoint:endpoint));  
sumV3 = sum(ww(startpoint:endpoint).*actualI(startpoint:endpoint));  %for only one R
sumV3ch = sum(ww(startpoint:endpoint).*Ich(startpoint:endpoint));  %for charge
sumV4di = sum(ww(startpoint:endpoint).*Idi(startpoint:endpoint));   %for discharge


sum1 =  sum(ww(startpoint:endpoint));

sumVacV1 = sum(ww(startpoint:endpoint).*actualV(startpoint:endpoint).*getV1(startpoint:endpoint));
sumVacV2 = sum(ww(startpoint:endpoint).*actualV(startpoint:endpoint).*getV2(startpoint:endpoint));
sumVacV3 = sum(ww(startpoint:endpoint).*actualV(startpoint:endpoint).*actualI(startpoint:endpoint));   %for only one R
sumVacV3ch = sum(ww(startpoint:endpoint).*actualV(startpoint:endpoint).*Ich(startpoint:endpoint));  %for charge
sumVacV4di = sum(ww(startpoint:endpoint).*actualV(startpoint:endpoint).*Idi(startpoint:endpoint));   %for discharge


sumVac = sum(ww(startpoint:endpoint).*actualV(startpoint:endpoint));


%Original before weight function

%%Notation is V1 for first alpha sequence, V2 for second, V3 for I sequence
%sumV1sq = sum(getV1(startpoint:endpoint).^2);
%sumV2sq=sum(getV2(startpoint:endpoint).^2);    %for second alpha
%sumV3sq = sum(actualI(startpoint:endpoint).^2);

%sumV1V2 = sum(getV1(startpoint:endpoint).*getV2(startpoint:endpoint));
%sumV1V3 = sum(getV1(startpoint:endpoint).*actualI(startpoint:endpoint));
%sumV2V3 = sum(getV2(startpoint:endpoint).*actualI(startpoint:endpoint));


%sumV1 = sum(getV1(startpoint:endpoint));
%sumV2 = sum(getV2(startpoint:endpoint));
%sumV3 = sum(actualI(startpoint:endpoint));


%sum1 =  length(getV1(startpoint:endpoint));  %N, the number of indicies

%sumVacV1 = sum(actualV(startpoint:endpoint).*getV1(startpoint:endpoint));
%sumVacV2 = sum(actualV(startpoint:endpoint).*getV2(startpoint:endpoint));
%sumVacV3 = sum(actualV(startpoint:endpoint).*actualI(startpoint:endpoint));
%sumVac = sum(actualV(startpoint:endpoint));

% What size of fit do we do? Two CPEs and an R or just one CPE



if (one_or_two_CPE==2)
    %two CPEs, one R
    M=[ sumV1sq sumV1V2  sumV1V3 sumV1 ;   sumV1V2 sumV2sq sumV2V3 sumV2;  ...
    sumV1V3 sumV2V3 sumV3sq sumV3;  sumV1 sumV2 sumV3 sum1 ] %matrix
    rhs = [sumVacV1 ; sumVacV2; sumVacV3;  sumVac]   %vector

    %Add small value to diagonal; hopefully helps with conditioning
    %small_val=200.0;
    %M=M+small_val*eye(4);   
    
    
    %M output = rhs;   %matrix equation
    output= inv(M)*rhs;  
     %output= pinv(M)*rhs;    %try dealing with near singluar cases

%   %12 dec 2022
%    %Fix to Vance's suggested.
%    disp('****FIXING VALUES*****')
%    output(1)=1/7500;
%    output(2)=1/50;
%    output(3)=0.15;
%    output(4)=4.000;

    CFfitted1 = 1/output(1);  
    CFfitted2 = 1/output(2); 
    Rfittedch = output(3);
    Rfitteddi = output(3);   %same
    constV=output(4);  

end
if (one_or_two_CPE==1)
    %one CPE, one R
   
    %if we do just one CPE ignore the second set. Ignore anything with a
    %V2...
    M=[ sumV1sq  sumV1V3 sumV1 ;    ...
    sumV1V3 sumV3sq sumV3 ; ...
    sumV1 sumV3 sum1];
    
    rhs=[sumVacV1; sumVacV3;  sumVac];   %vector

    %M output = rhs;   %matrix equation
    output= inv(M)*rhs;  
    
    CFfitted1 = 1/output(1);  
    CFfitted2 = inf; 
    Rfittedch = output(2); 
    Rfitteddi = output(2);   %same for discharge
    constV=output(3);  
    
    disp('***Fitting only the first CPE***')
end

if (one_or_two_CPE==3)
    %two CPEs, R charge and R discharge
    M=[ sumV1sq   sumV1V2   sumV1V3ch   sumV1V4di   sumV1 ;  ...
        sumV1V2   sumV2sq   sumV2V3ch   sumV2V4di   sumV2;  ...
        sumV1V3ch sumV2V3ch sumV3chsq   sumV3chV4di sumV3ch;  ...
        sumV1V4di sumV2V4di sumV3chV4di sumV4disq   sumV4di; ...
        sumV1     sumV2     sumV3ch     sumV4di     sum1 ]      %matrix
    rhs = [sumVacV1 ; sumVacV2; sumVacV3ch;  sumVacV4di;  sumVac]   %vector

    %Add small value to diagonal; hopefully helps with conditioning
    %small_val=200.0;
    %M=M+small_val*eye(4);   
    
    
    %M output = rhs;   %matrix equation
    output= inv(M)*rhs;  
     %output= pinv(M)*rhs;    %try dealing with near singluar cases

%   %12 dec 2022
%    %Fix to Vance's suggested.
%    disp('****FIXING VALUES*****')
%    output(1)=1/7500;
%    output(2)=1/50;
%    output(3)=0.15;
%    output(4)=4.000;

    CFfitted1 = 1/output(1);  
    CFfitted2 = 1/output(2); 
    Rfittedch = output(3);
    Rfitteddi = output(4);   %same
    constV=output(5);  

end



 

if (~fit_const)
    %if we don't fit the constant voltage, work out fit with a 3x3 matrix
    M = [M(1:3,1:3)];
    rhs = rhs(1:3)-actualV(1)*[sumV1; sumV2; sumV3];
    output=inv(M)*rhs;
    
     %12 dec 2022
%    %Fix to Vance's suggested.
%    disp('****FIXING VALUES*****')
%    output(1)=1/7500;
%    output(2)=1/50;
%    output(3)=0.15;
%    output(4)=4.000;
    
     CFfitted1 = 1/output(1);  
    CFfitted2 = 1/output(2); 
    Rfittedch = output(3);    
    Rfitteddi = output(3); %same
   constV=actualV(1);    %constant is what we specify it to start
   
   
  
   
   
end

if (~suppress)
 disp(['CF1 is ' num2str(CFfitted1)]);
 disp(['CF2 is ' num2str(CFfitted2)]);
 disp(['Rs_ch is ' num2str(Rfittedch)]);
 disp(['Rs_di is ' num2str(Rfitteddi)])
 disp(['Const V is ' num2str(constV)]);
end
 
%Vout = output(1)*getV1 + output(2)*getV2 + output(3)*actualI + constV;    %the fitted output vector
Vout = getV1/CFfitted1 + getV2/CFfitted2 + Rfittedch*Ich + Rfitteddi*Idi + constV;   %better to lose the reference to output index - more versatile.

%17 January 2023
%Mask out anything below 3.5 V in the calculation

accept_points1=(Vout(startpoint:endpoint)>3.5);  %identify which points are okay
accept_points2=(Vout(startpoint:endpoint)<4.1);
accept_points=accept_points1.*accept_points2;  %Gendrive should have all accepted within 3.5-4.1 V

difference_sequence = Vout(startpoint:endpoint) - actualV(startpoint:endpoint);
S = sum((difference_sequence.*accept_points).^2);

lsq=S/sum(accept_points);  %value per accepted point




%old case
%S=sum(  (  Vout(startpoint:endpoint) - actualV(startpoint:endpoint)  ).^2  );
%lsq = S/length(actualV(startpoint:endpoint)) ;   %value per point

% %Try reporting the gradient diff
% diffVout=(Vout(2:end) - Vout(1:end-1));
% diffactualV=(actualV(2:end) - actualV(1:end-1));
% lsq=sum( (diffVout - diffactualV).^2 )/length(diffVout); 

if (~suppress)   %show the output
 disp(['lsq is ' num2str(lsq)])

 delta_a1 = sqrt(S/sumV1sq);   %error in 1/CF1
 deltaCF1 = CFfitted1^2*delta_a1;

 delta_a2 = sqrt(S/sumV2sq);   %error in 1/CF2
 deltaCF2 = CFfitted2^2*delta_a2;

 delta_b = sqrt(S/sumV3sq);  %error in Rch
 deltaRs = delta_b;

 delta_c = sqrt(S/(endpoint-startpoint));  %error in constant voltage

 %error in fit (p8050)
 disp(['error in CF1 is ' num2str(deltaCF1) ])
 disp(['error in CF2 is ' num2str(deltaCF2) ])
 disp(['error in Rs is ' num2str(deltaRs) ])
 disp(['error in const is ' num2str(delta_c) ])
end

end