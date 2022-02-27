% Supplementary material for the paper:
% Granata et al: " Integration of transcriptomic data in a genome-scale metabolic 
% model to investigate the relation between obesity and breast cancer"

% this code makes use of the Lee-12 matlab code, provided on request by the authors 
% see Lee et al.,  "Improving metabolic flux predictions using absolute gene expression data"

% set the solvers
changeCobraSolver('gurobi6', 'LP');
changeCobraSolver('gurobi6', 'MILP');
% path of the SBML model name
model_name      = 'iAdipocytes1809.xml';
% comment/uncomment the desired expression data
% lean
genedata_name = 'lean.txt'; 
% obese
%genedata_name = 'obese.txt'; 

%%%%%% now predict flux distribution (v_gene_exp) based on the provided expression data

% load model
model       = importModel(model_name);

% load transcript data
genedata	= importdata(genedata_name);

if size(genedata.data,2) == 3
    genenames	= strtrim(cellstr(num2str(genedata.data(:,1))));
    gene_exp	= genedata.data(:,2);
    gene_exp_sd	= genedata.data(:,3);
else
    genenames	= genedata.textdata(:,1);
    genenames(1)= [];
    gene_exp	= genedata.data(:,1);
    gene_exp_sd	= genedata.data(:,2);    
end

% map gene weighting to reaction weighting
disp('Mapping gene to reactions...');
[rxn_exp,rxn_exp_sd] = geneToReaction(model,genenames,gene_exp,gene_exp_sd);
disp('...Done!')
rxn_exp_sd(rxn_exp_sd == 0) = min(rxn_exp_sd(rxn_exp_sd>0))/2;

% We extract the specific metabolic model by setting up the upper and lower bounds of TAG and Glucose uptake fluxes 
% in the two specific conditions.

% Here we set up the input fluxes using  experimentally measured values for TAG and Glucose uptake (see the paper for details).
% comment/uncomment to choose lean/obese constraints, according to the data imported above. 
% Choose the time point you want (either tp=4 or tp=5) 

% "HMR_n" is the reaction number corresponding to the specific flux we want to constrain
TAG_flux=find(strcmp('HMR_9023',model.rxns)); 
glucose_uptake=find(strcmp('HMR_9034',model.rxns));

% TAG EXTRACTION:
% Time point tp=4
%lean upper and lower bound for TAG 
model.lb(TAG_flux)	= 0.345519;	
model.ub(TAG_flux)	= 0.345519;
%obese upper and lower bound for TAG
%model.lb(TAG_flux)	= 0.1028331;
%model.ub(TAG_flux)	= 0.1028331;

% Time point 5
%lean upper and lower bound for TAG
%model.lb(TAG_flux)	= 0.2481715;	
%model.ub(TAG_flux)	= 0.2481715;
%obese upper and lower bound for TAG
%model.lb(TAG_flux)	= 0.1546842;
%model.ub(TAG_flux)	= 0.1546842;

% GLUCOSE UPTAKE:
% Time point tp=4
%lean upper and lower bound for TAG 
model.lb(glucose_uptake)	= 1.585;
model.ub(glucose_uptake)	= 1.585;
%obese upper and lower bound for TAG 
%model.lb(glucose_uptake)	= 0.5715;
%model.ub(glucose_uptake)	= 0.5715;

% Time point tp=5
%lean upper and lower bound for TAG 
%model.lb(glucose_uptake)	= 1.1263;
%model.ub(glucose_uptake)	= 1.1263;
%obese upper and lower bound for TAG 
%model.lb(glucose_uptake)	= 0.7464;
%model.ub(glucose_uptake)	= 0.7464;

% Now run the gene expression constraint FBA, and get the resulting flux distribution
v_gene_exp      = dataToFlux(model,rxn_exp,rxn_exp_sd); 

