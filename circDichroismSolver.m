%%File name: CircularDichroism.m
%%Author: Quintin T. Davis
%%Data extraction, verification, quality control: Lauren R. Johnson
%%completed as part of MOLB 4180, and submitted as part of a report on circular
%%dichroism for BE 4820, University of Wyoming
%%Description: Data come from MOLB 4180 proteins laboratory course

%%Version 1.0 simply added weighted components and user judgement
%%ascertained best fit
%%Version 1.1 added a calculation for mean square error, and optimized 
%%  by reiteration for a minimal MSE error
%%Version 1.3 added an algorithm to automatically extrapolate the data
%%for reference spectra, and commenting was improved.
%Revised 4/24/2010
%Revised 4/30/2010
%Revised 5/2/2012
clear all;
%creates a vector for x-axis, wavelength of light (nm)
lemp = [190,191,192.5,195,197,200,202,205,208,210,211,214,215,217,220,222,225,230,234,238,240,250]; %Greenfield and Fasman Data
lones = (190:1:250); %to match resolution (1 nm increments) from class data
%%Secondary structure standards created by visual observation of handout and normalized.%%
bShet1 = [22400,25300,30000,31900,30000,24300,19300,5700,-4700,-10800,-12100,-16400,-17900,-18400,-15700,-13800,-11400,-6400,3600,-1400,700,0];
aHelix1 = [74800,76900,73300,64300,44300,14300,0,-25000,-32600,-32400,-32100,-31000,-31400,-33100,-35300,-35700,-32400,-21900,-11400,-4300,-3300,0];
coil1 = [-32200,-34700,-37500,-41000,-41900,-36400,-25600,-14500,-3400,-1400,0,3500,4100,4600,4400,3900,2700,800,0,-140,-150,0];
%above data from Fasman and Greenfield 1969

%counters for empirical index (k) and by ones (c)
c = 1;
k = 1;
%this loop extrapolates between data obtained from Greenfield et al [2].to generate full dataset
while c < length(lones)
    if lones(c) == lemp(k)
        aHelix(c) = aHelix1(k);
        bShet(c) = bShet1(k);
        coil(c) = coil1(k);
        c = c+1;
    else
        while lones(c) < lemp(k)
               aHelix(c) = (((lemp(k)-lemp(k-1)) - (lemp(k)-lones(c)))/ (lemp(k)-lemp(k-1))) * (aHelix1(k)-aHelix1(k-1)) + aHelix1(k-1);
               bShet(c) = (((lemp(k)-lemp(k-1)) - (lemp(k)-lones(c)))/ (lemp(k)-lemp(k-1))) * (bShet1(k)-bShet1(k-1)) + bShet1(k-1);
               coil(c) = (((lemp(k)-lemp(k-1)) - (lemp(k)-lones(c)))/ (lemp(k)-lemp(k-1))) * (coil1(k)-coil1(k-1)) + coil1(k-1);
               c = c+1;
        end
    end
    k = k+1;
end
aHelix(c) = 0;
bShet(c) = 0;
coil(c) = 0;
alpha = aHelix / max(abs(aHelix));
beta = bShet / max(abs(aHelix)); %normalize all values to the greatestmagnitude of alpha Helix signal
coill = coil / max(abs(aHelix));
%%The BSA sample data was constructed from a graph (no numerical data was given by instructor) by visual
%%observation and also normalized.The spectra from all 3 sample conditions can be solved for by uncommenting the 
%%desired 'BSASpecTest' variable below.
%BSASpecTest = [10 13.75 17.5 18.5 19.5 19.75 20.0 17.5 15 10.5 6 1 -4 -6.5 -9.75 -12.5 -15 -16.25 -17.5 -17.75 -18.0 -17.5 -17.0 -16.5 -16.0 -15.9 -15.8 -15.4 -15.0 -14.5 -14.0 -13.5 -13.0 -12.5 -12 -10 -8 -7 -6.3 -5.1 -4.0 -3.1 -2.2 -1.9 -1.6 -1.3 -1 -0.9 -0.8 -0.65 -0.5 -.375 -0.25 -.18 -0.13 -0.09 -.05 0 0 0 0];
%the above is for the BSA in TFE
%BSASpecTest = [10 11.75 13.5 14.5 15.5 15.25 15 14.25 12.5 10.75 9 7.25 6.5 3.25 0 -2.5 -5 -7 -9 -9.7 -10.4 -10.375 -10.35 -10.2 -10.05 -10.025 -10 -9.6 -9.2 -9.1 -9 -8.6 -8.2 -7.6 -7 -6.5 -6 -5 -4.0 -3.875 -3.75 -2.9 -2.05 -2.025 -2 -1.9 -1.8 -1.775 -1.75 -1.74 -1.73 -1.715 -1.7 -1.7 -1.7 -1.7 -1.7 -1.7 0 0 0];
%the above is the hot BSA spectrum (80 degrees)
BSASpecTest = [14.5 15.25 16 16.5 17 17 17 16.5 16 13 10 4 -2 -4 -6 -8.5 -11 -11.25 -11.5 -11.625 -11.75 -11.625 -11.5 -11.25 -11 -10.75 -10.5 -10.375 -10.25 -10 -9.75 -9.5 -8.75 -8 -7 -6 -5.25 -4.5 -4 -3.5 -2.75 -2 -1.75 -1.5 -1 -.5 -.5 -.5 -.5 -.5 -.5 -.5 -.5 -.5 -.5 -.5 -.5 0 0 0 0 ];
%above data is from BSA w/o TFE @ room temp.
BSASpec = BSASpecTest / max(BSASpecTest); % normalize
user = 'y'; %y/n to continue is used later
tolerance = input('Enter the highest desired mean square error');
iter = input('Enter the number of algorithm iterations desired');
SqE = (10000); %throwaway SqE set unreasonably high to start
%Weights are equal to start with at 33%. These correspond to each secondary structures contribution to the CD spectrum
awtemp = 1/3;
bwtemp = 1/3;
cwtemp = 1/3;
i = 0; %iteration counter
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock))); %initialize Mersenne twister pseudo-random number gen.
while user == 'y' & i < iter
%%random weight generated for each type of secondary structure. total
%%weight stored to determine percentage later. 20% randomized each time, best weights retained.%%
    awt = 0.8 * awtemp + .2 * rand ;
    bwt = 0.8 * bwtemp + .2 * rand ;
    cwt = 0.8 * cwtemp + .2 * rand ;
    Test = awt * alpha + bwt * beta + cwt * coill;
    Test = Test / max(abs(Test));
    k = 1; %counter for computing square error
    SqETemp = 0; %This will be used as an accumulator, so is set to zero to begin.
    while k < 62
        SqETemp = SqETemp + (BSASpec(k) - Test(k))^2;
        k = k + 1;
    end
    if SqETemp < SqE,
        SqE = SqETemp;
        awtemp = awt;
        bwtemp = bwt;
        cwtemp = cwt;
    end
    MSqE = SqE/length(lones); %divide by number of datapoints to get average
    if MSqE < tolerance %ask user if they want to continue once threshold MSE is reached.
        disp('current best mean square error is'), disp(MSqE)
        disp('current tolerance is '), disp(tolerance)
        user = input('continue y/n?','s');
        if user == 'y'
            tolerance = input('Enter new value for MSE tolerance');
        end
    end

    i = i + 1; %increment iteration counter
end
%Display the results of the algorithm as % contribution of each type of structure, then plot the result for visual
%verification. 
tot = awtemp + bwtemp + cwtemp;
Test = awtemp * alpha + bwtemp * beta + cwtemp * coill;
Test = Test / max(abs(Test));
plot(lones,BSASpec,'b-',lones,BSASpec,'go',lones,Test,'r-',lones,Test,'gx');
xlabel('wavelength (nm)');
ylabel('normalized magnitude');
title('BSA Spectrum (o) and simulated CD results (x)')
sprintf('The best match for secondary structure contributions to the sample CD spectra, achieved in ',iter,' iterations or within user defined mean square error tolerance is ')
disp('alpha helix contribution %:');disp(100*awtemp/tot);
disp('beta sheet contribution %:'); disp(100*bwtemp/tot);
disp('random coil contribution %:');disp(100*cwtemp/tot);
disp('the mean square error for the result of this composition is');disp(MSqE)
