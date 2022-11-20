%% WORKSPACE
cd 'C:\github\Panyawarin';
initCobraToolbox;
%setRavenSolver('cobra');
%% Load the model
load('ComplementaryData/iPS1469.mat');
Modelori = model
%% STEP4 : MODEL VALIDATION
% The model validation was performed against dark experimentation by Raethong et al. (2018) 
% and light-programming experimentation by Thananusak et al. (2020).

    %%% Dark condition 
    % The growth simulation on glucose cultures grown under dark conditions.
    
iPC1469Raven_dark = Modelori;
iPC1469Raven_dark = setParam(iPC1469Raven_dark,'lb',{'bmOUT','cordycepinOUT'},[0 0]);
iPC1469Raven_dark = setParam(iPC1469Raven_dark,'eq',{'matp'},1);
iPC1469Raven_dark = setParam(iPC1469Raven_dark,'ub',{'bmOUT','cordycepinOUT'},[1000 1000]);
iPC1469Raven_dark = setParam(iPC1469Raven_dark,'obj',{'bmOUT'},1);
iPC1469Raven_dark = setParam(iPC1469Raven_dark,'eq',{
        'neurosporaxanthin','rt_R0001','rt_R0002','rt_R0003','rt_R0004','rt_R0005','rt_R0007','R09782_c','rt_R0009'},[
           0 0 0 0 0 0 0 0 0]);%block caroteniods ractions
sold = solveLP(iPC1469Raven_dark,1);

sugars = {'glcIN''sucIN'};
uptake = [0.1448,0.0815]; %observed from experiment (Raethong  et al., 2020)   
sumax(11) = sold;
printFluxes(iPC1469Raven_dark, sumax(11).x,true);

for i = 1:numel(uptake)
    model=setParam(iPC1469Raven_dark,'ub',sugars,[0 0]);
    model=setParam(model,'ub',sugars(i),uptake(i));
    sugar = sugars(i);
    sold = solveLP(model,1);
    sumax(i) = sold;
    umax = num2str(sold.f*-1);
    fprintf([char(sugar) '\t' num2str(sold.f*-24) ' per day' '\n']);
end

finalValidateModel0 = setParam(iPC1469Raven_dark,'obj',{'cordycepinOUT'},1);
finalValidateModel1 = setParam(iPC1469Raven_dark,'obj',{'exc sphinganine'},1);
finalValidateModel3 = setParam(iPC1469Raven_dark,'obj',{'exc phytosphingosine'},1);
finalValidateModel4 = setParam(iPC1469Raven_dark,'obj',{'exc sphingosine'},1);
finalValidateModel5 = setParam(iPC1469Raven_dark,'obj',{'exc neurosporaxanthin'},1);

sol0 = solveLP(finalValidateModel0,1);
fprintf(['production rate of cordycepin' '\t' num2str(sol0.f*-1) ' mmol/g DW/h' '\n']);
sol1 = solveLP(finalValidateModel1,1);
fprintf(['production rate of sphinganine' '\t' num2str(sol1.f*-1) ' mmol/g DW/h' '\n']);
sol3 = solveLP(finalValidateModel3,1);
fprintf(['production rate of phytosphingosine' '\t' num2str(sol3.f*-1) ' mmol/g DW/h' '\n']);
sol4 = solveLP(finalValidateModel4,1);
fprintf(['production rate of sphingosine' '\t' num2str(sol4.f*-1) ' mmol/g DW/h' '\n']);
sol5 = solveLP(finalValidateModel5,1);
fprintf(['production rate of neurosporaxanthin' '\t' num2str(sol5.f*-1) ' mmol/g DW/h' '\n']);

    %%% light-programming condition
     % The growth simulation on glucose cultures grown under Light condition
    
iPC1469Raven_light = Modelori   
iPC1469Raven_light = setParam(iPC1469Raven_light,'obj','bm_light',1); 
iPC1469Raven_light = setParam(iPC1469Raven_light,'lb',{'bm_light','cordycepinOUT'},[0 0]);
iPC1469Raven_light = setParam(iPC1469Raven_light,'ub',{'bm_light','cordycepinOUT'},[1000 1000]);
iPC1469Raven_light = setParam(iPC1469Raven_light,'eq',{'bmOUT'},0);
solL = solveLP(iPC1469Raven_light,1);

sugars = {'glcIN''sucIN'};
uptake = [0.0116,0.0122]; %observed from experiment(Thananusak et al., 2020)
sumax(1) = solL;
printFluxes(iPC1469Raven_light, sumax(1).x,true);

for i = 1:numel(uptake)
    model=setParam(iPC1469Raven_dark,'ub',sugars,[0 0]);
    model=setParam(model,'ub',sugars(i),uptake(i));
    sugar = sugars(i);
    soll = solveLP(model,1);
    sumax(i) = sold;
    umax = num2str(sold.f*-1);
    fprintf([char(sugar) '\t' num2str(soll.f*-24) ' per day' '\n']);
end


finalValidateModel6 = setParam(iPC1469Raven_light,'obj',{'cordycepinOUT'},1);
finalValidateModel7 = setParam(iPC1469Raven_light,'obj',{'exc sphinganine'},1);
finalValidateModel8 = setParam(iPC1469Raven_light,'obj',{'exc phytosphingosine'},1);
finalValidateModel9 = setParam(iPC1469Raven_light,'obj',{'exc sphingosine'},1);
finalValidateModel10 = setParam(iPC1469Raven_light,'obj',{'exc neurosporaxanthin'},1);

sol6 = solveLP(finalValidateModel6,1);
fprintf(['production rate of cordycepin' '\t' num2str(sol0.f*-1) ' mmol/g DW/h' '\n']);
sol7 = solveLP(finalValidateModel7,1);
fprintf(['production rate of sphinganine' '\t' num2str(sol1.f*-1) ' mmol/g DW/h' '\n']);
sol8 = solveLP(finalValidateModel8,1);
fprintf(['production rate of phytosphingosine' '\t' num2str(sol3.f*-1) ' mmol/g DW/h' '\n']);
sol9 = solveLP(finalValidateModel9,1);
fprintf(['production rate of sphingosine' '\t' num2str(sol4.f*-1) ' mmol/g DW/h' '\n']);
sol10 = solveLP(finalValidateModel10,1);
fprintf(['production rate of neurosporaxanthin' '\t' num2str(sol5.f*-1) ' mmol/g DW/h' '\n']);
    
    
 %%%%%%%%%%%%%%%%%%%%%%%%

 %                                         Uptake rate      |       Growth rate,μmax(h-1)      |  Error rate %
 %                                        (mmol gDW-1 h-1)  |   Experiments     |   Prediction |  
 %                                                          |                   |              |         
%Light-programming condition of Glucose (a) 0.1961 ± 0.0568 |   0.0116 ± 0.0040 |  0.0116      |   0.26                   
%Light-programming condition of Sucrose (a) 0.1428 ± 0.0349 |   0.0122 ± 0.0017 |              |                                  
%Dark condition of Glucose (b)              0.1448 ± 0.0872 |   0.0100 ± 0.0027 |  0.0100      |   0.40                                             
%Dark condition of Sucrose (b)              0.0815 ± 0.0250 |   0.0114 ± 0.0014 |  0.0117      |   2.42            
                                                        
% (a) Thananusak et al., 2020
% (b) Raethong  et al., 2020   
   
  
 %                                 %Carotenoid content (mg g DCW−1) %
%Light-programming condition of Glucose (a) 1.4692 ± 0.0122                       
%Light-programming condition of Sucrose (a) 1.4590 ± 0.0052                                   
%Dark condition of Glucose (a)              0.0244 ± 0.0051                                                
%Dark condition of Sucrose (a)              0.0234 ± 0.0095       

% (a) Thananusak et al., 2020    
    
    
    
    
    
 
    
    
    
    
    
    
    
    
    
    
    
    