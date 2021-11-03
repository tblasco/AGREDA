# AGREDA
AGREDA (AGORA-based REconstruction for Diet Analysis) is a new repository of genome-scale metabolic models of the human gut microbiota. AGREDA includes degradation pathways for many diet derived compounds, which are mainly metabolized by the gut microbiota. In particular, AGREDA incorporates 179 degradation pathways of phenolic compounds, which play an important role in human health and nutrition, and are closely related to the human gut microbiota.

AGREDA mixed-bag model and subsequent species models can be found at **AGREDA_v1.0.0.zip** file.

For further information, please refer to:
* Telmo Blasco	tblasco@tecnun.es
* Francesco Balzerani	fbalzerani@tecnun.es
* Iñigo Apaolaza	iaemparanza@tecnun.es
* Francisco J. Planes	fplanes@tecnun.es

## Citing AGREDA
Blasco, T., Pérez-Burillo, S., Balzerani, F., Hinojosa-Nogueira, D., Lerma-Aguilera, A., Pastoriza, S., ... & Planes, F. J. (2021). An extended reconstruction of human gut microbiota metabolism of dietary compounds. Nature communications, 12(1), 1-12.

# PIPELINE

Please refer to the following sections for a pipeline execution description.

## SOFTWARE REQUIREMENTS

* It will be necessary to install Matlab. We encourage users to download release R2018a.

* Free academic licenses for the IBM CPLEX solver can be obtained from https://www.ibm.com/developerworks/community/blogs/jfp/entry/CPLEX_Is_Free_For_Students?lang=en. We    encourage users to download CPLEX version 12.8.0.

* Cobratoolbox github link https://github.com/opencobra/cobratoolbox.

* In order to generate the supra-organism model, please download AGORA 1.03 version (with mucins) from https://vmh.life/#downloadview and locate metabolic reconstructions in .mat format at **./pipeline/input/AGORA_models/** directory.

* In order to run any section of the pipeline, please uncompress **AGREDA_UNBALANCED.zip**, **modelSEED.zip** and **supraOrganism-all_models.zip** files located at 
**./pipeline/output/Models/** directory.

## EXECUTION

Different sections can be run independently by **main.m** script as all the necessary data is located at input directory.

* **main.m Section 1:** AGORA mixed-bag model -> Building AGORA species level and mixed-bag models.
* **main.m Section 2:** AGREDA building -> Building AGREDA mixed-bag model.
* **main.m Section 3:** AGREDA balancing -> Single species analysis

## EXECUTION TIME 

* **main.m Section 1:** AGORA mixed-bag model -> about 1,5 hours.
* **main.m Section 2:** AGREDA building -> about 9 hours.
* **main.m Section 3:** AGREDA balancing -> about 2 hours

# FOLDER CONTENT

Please refer to the following sections for a pipeline execution description.

## INPUT FOLDER

#### AGORA_models:

Folder containing AGORA models in .mat format.

#### Annotation:

  * **AGORA_SEED_with_exc.xlsx:** link between AGORA and modelSEED metabolites.
  * **Input_metabolites.xlsx:** i-Diet metabolites/nutrients list.
  * **manual_added_ECnumbers_to_rxns.xlsx:** manually added EC numbers to modelSEED metabolic network reactions.
  * **manual_added_trRules_to_ECnumbers.xlsx:** manually added taxonomic information to EC numbers.
  * **manual_added_ECnumbers_to_rxns.xlsx:** manually added taxonomic information to modelSEED metabolic network reactions.
  * **mergeMetabolitesSEED.xlsx:** common metabolites from modelSEED that have to be merged.
  * **new_metabolites.xlsx:** metabolites added from expert knowledge.
  * **new_reactions.xlsx:** reactions added from expert knowledge.
  * **no_bacteria.xlsx:** reactions regarding species not present in the model.
  * **removeMetabolitesSEED.xlsx:** metabolites with limited evidence that have to be removed from the model.
  * **speciesInfo.xlsx:** information of species present in AGORA.
  * **table_similarity.xlsx:** repeated metabolites inside modelSEED.
  * **targetMetabolites.xlsx:** set of metabolites to add a transport reaction, based on the evidence found at HMDB.

##### EC_Numbers:

  * **myRAST:** Folder containing EC numbers obtained from myRAST.
  * **KEGG:** Folder containing EC numbers obtained from KEGG database.

##### FVA_Species: 

Folder containing fastFVA information during species balancing process.

##### OptionsFastcore:

  * **blockedDietRxnsOut_AGREDA_UNBALANCED.mat:** blocked i-Diet metabolites production reaction.
  * **blockedDietRxnsUp_AGREDA_UNBALANCED.mat:** blocked i-Diet metabolites uptake reaction.
  * **fastFVAAgoraSeedExpertBounded_AGREDA_UNBALANCED.mat:** fastFVA result of global merged model.
  * **removeNotAnnotatedRxnsFastFVA.mat:** fastFVA result after removing reactions without taxonomy.

##### SpeciesToMerge:
  
  * **all_models.mat:** All the species present in AGORA.

## CODE FOLDER

  * **addExchangesInfo.m:** function to add transport reaction to a set of metabolites based on the evidence found at HMDB.
  * **addMetabolites.m:** function to add the metabolites from expert knowledge.
  * **addReactions.m:** function to add the reactions from expert knowledge.
  * **applyFastcoreDiet.m:** function to integrate i-Diet metabolites efficiently to generate AGREDA model, based on fastCoreWeighted.
  * **balanceSpecies.m:** function to balance AGREDA model based on the species paths.
  * **ConcatenateModels.m:** function to concatenate AGORA and SEED metabolic networks with all their fields.
  * **Create_exchange_reaction_agora.m:** function to create exchange reactions (if necessary) for i-Diet metabolites present in AGORA.
  * **Create_exchange_reaction_seed.m:** function to create exchange reactions (if necessary) for i-Diet metabolites present in SEED.
  * **FBA.m:** function to apply Flux Balance Analysis.
  * **manageEC_pipeline.m:** function to add taxonomic information to SEED metabolic network through EC number information.
  * **mapToAGREDA.m:** function to contextualize species information to the AGREDA model.
  * **Merge_mets_by_structure.m:** function to merge information of repeated metabolites in SEED.
  * **mergeDupRxns:** function to merge duplicated reactions in a model.
  * **mergeMets.m:** function to merge metabolite information from modelSEED.
  * **mergeSEEDmetabolites.m:** function to merge AGORA and SEED models.
  * **mergeSpeciesLevel.m:** function to concatenate all AGORA reconstructions.
  * **mergeSupraOrganism.m:** function to merge the species-level model into a mixed-bag model.
  * **removeNotAnnotatedRxns.m:** function to remove not annotated rxns of the final model and apply a model reconstruction with fastFVA.
  * **removeProblematicMets.m:** function to remove metabolites from modelSEED with limited evidence.
  * **removeRxnsAndMets.m:** function to remove a set of reactions and their related metabolites (if necessary) from a model, based on a reaction index.

## OUTPUT FOLDER

##### Models:

Folder containing the models generated during the different steps of the pipeline.
