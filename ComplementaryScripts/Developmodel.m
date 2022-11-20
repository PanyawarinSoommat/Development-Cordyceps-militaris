%% Development of regulatory genome-scale modeling of Cordyceps militaris for probing light-responsive mechanism
% Purposed:
% The process was carried through the following steps:
% STEP 1: TEMPLATE MODEL PREPARATION
% STEP 2: ADD NEWLY-IDENTIFIED REACTIONS
% STEP 3: ADD BIOMASS REACTIONS
% STEP 4: MODEL VALIDATION
%% WORKSPACE
cd 'C:\github\Panyawarin'
% Initialize COBRA toolbox
initCobraToolbox;
%setRavenSolver('cobra');
%% STEP 1: TEMPLATE MODEL PREPARATION
% DVELOPMENT OF THE MODEL
 %  The Genome-scale metabolic models (GSMMs) of C. militaris, iPC1469 was used as a template to improved genome-scale
 % metabolic network reconstructions.
 
model=importModel('ComplementaryData\revised_iPC1469.xml') % Import the model

iPC1469Model = model;
iPC1469Cobra = ravenCobraWrapper(iPC1469Model);
iPC1469Raven = ravenCobraWrapper(iPC1469Cobra);

%Check the iPC1469 model 
finalValidateModel = iPC1469Raven;
finalValidateModel = setParam(finalValidateModel,'lb',{'bmOUT','cordycepinOUT'},[0 0]);
finalValidateModel = setParam(finalValidateModel,'eq',{'matp'},1);
finalValidateModel = setParam(finalValidateModel,'ub',{'bmOUT','cordycepinOUT'},[1000 1000]);
finalValidateModel = setParam(finalValidateModel,'obj',{'bmOUT'},1);
sol1 = solveLP(finalValidateModel,1);

sugars = {'glcIN' 'fruIN' 'arabIN' 'xylIN' 'sucIN'};
uptake = [0.1448,0.1472,0.1074,0.0681,0.0815]; %observed from experiment
sumax(11) = sol1;
printFluxes(finalValidateModel, sumax(11).x,true);


for i = 1:numel(uptake)
    model=setParam(finalValidateModel,'ub',sugars,[0 0 0 0 0]);
    model=setParam(model,'ub',sugars(i),uptake(i));
    sugar = sugars(i);
    sol = solveLP(model,1);
    sumax(i) = sol;
    umax = num2str(sol.f*-1);
    fprintf([char(sugar) '\t' num2str(sol.f*-24) ' per day' '\n']);
end

%% STEP 2: ADD NEWLY-IDENTIFIED REACTIONS
% 2.1 Addition of new metabolites to the templates network 
 % New metabolites listed in NEW_Metabolites were introduced to the templates network by addMets function.

[~, newMets]=xlsread('ComplementaryData/supplementary_new.xlsx','NEW_Metabolites');
metsToAdd = struct();
metsToAdd.mets = newMets(2:end,1);
metsToAdd.metNames = newMets(2:end,2);
metsToAdd.metFormulas = newMets(2:end,3);
metsToAdd.compartments = 'c';
metsToAdd.inchis = newMets(2:end,9);
addMetsModel=addMets(iPC1469Raven,metsToAdd);
addMetsModel = contractModel(addMetsModel);

addMetsModel.equations = constructEquations(addMetsModel);
solmet = solveLP(addMetsModel,1);
solmet.f

%exportForGit(addMetsModel,'modelmet','',{'xlsx'});

% 2.2 Addition of new reactions to the templates network 
[~, SheetS]=xlsread('ComplementaryData/supplementary_new.xlsx','NEW_Reactions');
newRxns = struct();
newRxns.rxns = SheetS(2:end,1);
newRxns.rxnNames = SheetS(2:end,2);
newRxns.equations = SheetS(2:end,3);
newRxns.eccodes = SheetS(2:end,4);
newRxns.grRules = SheetS(2:end,5);
newRxns.rxnNotes = SheetS(2:end,13);

addRxnModel = addRxns(addMetsModel,newRxns,3,'',true,true);
addRxnModel.equations = constructEquations(addRxnModel);
addRxnModel = sortModel(addRxnModel);
model_new1 = contractModel(addRxnModel);

model_new1.lb(isinf(-addRxnModel.lb))=-1000; 
model_new1.ub(isinf(addRxnModel.ub))=1000;

model_new1.equations = constructEquations(model_new1);
solrxn = solveLP(model_new1,1);
printFluxes(model_new1, solrxn.x,true);
iPC1469Raven1=model_new1;


%% STEP3 : ADD BIOMASS REACTIONS
% The template model does not contain light biomass reactions so biomass reactions must be added related to the light.
[~, SheetS]=xlsread('ComplementaryData/supplementary_new.xlsx','NEW_Biomass');
add_Biomass = struct();
add_Biomass.rxns = SheetS(2:end,1);
add_Biomass.rxnNames = SheetS(2:end,2);
add_Biomass.equations = SheetS(2:end,3);
add_Biomass.rxnNotes = SheetS(2:end,13);

model_bm = addRxns(iPC1469Raven1,add_Biomass,3,'',true,false);
model_bm.equations = constructEquations(model_bm);
model_bm = sortModel(model_bm);
model_bm1 = contractModel(model_bm);


[a, b]=ismember(model_bm1.rxns,add_Biomass.rxns);
I=find(a);
model_bm1.rev(I)=0;
reducedModel1=model_bm1;
for i = 1:numel(add_Biomass.rxns)
    reducedModel1=setParam(reducedModel1,'lb',add_Biomass.rxns(i),[0]);
    reducedModel1=setParam(reducedModel1,'ub',add_Biomass.rxns(i),[1000]);
end

reducedModel1.equations = constructEquations(reducedModel1);
solrxn = solveLP(reducedModel1,1);
printFluxes(reducedModel1, solrxn.x,true);
iPC1469Raven2=reducedModel1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iPS1469 = iPC1469Raven2;
iPS1469.id = 'iPS1469';
iPS1469.name = 'iPS1469';
iPS1469.description = 'Development C. militaris GSMM ';
iPS1469.annotation = [];
printModelStats(iPS1469);

%%%%%%%%%%%%%%%%
%exportForGit(iPS1469,'iPS1469','',{'mat', 'txt', 'xlsx', 'xml'});


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
    