### Data and code for JEB paper on evolution of habitat and sexual isolation ######
### Title: The role of divergent host use and geography in the evolution of habitat isolation and sexual isolation among sister species of Belonocnema gall wasps #####
###### Authors: Linyi Zhang, Glen Ray Hood, James R. Ott, and Scott P. Egan ###
###### For questions about data and code, please contact Linyi Zhang at linyi.zhang@gwu.edu or linyizhangecnu@gmail.com 

**This repository includes R code file "data analysis code.R" for analyzing datasets "Bt host preference new.csv" and "mate.preference.csv"**

**The R code file "data analysis code.R" includes all the packages you need to load and all the steps you need to regenerate the data analysis and figures.**

+ **"Bt host preference new.csv" is the data file document the measurements of host preference across all three Belonocnema species.**

  - Column "Method" stands for the method of conducting host preference assay (cup or petri dish)."HP" stands for which host plant the wasp individual come from. "Qv-F" means waps invidividuals are from Q. vriginiana plants in sympatry. "Qv-T" means waps invidividuals are from Q. vriginiana plants in allopatry. "Qg" means waps invidividuals are from Q. geminata plants. 

  - Columns labeled from "0" to "20" represent the time interval the observation took. For instance, "0" represents observations of behavior before timing starts. "1" represents observation of behavior at the first time interval. 

+ **"mate.preference.csv" is the data file document the measurements of mate preference across all three Belonocnema species.**
- The column labeled "Method" indicates the specific method used to conduct the host preference assay, which can be either a cup or a petri dish. The abbreviation "HP" is used to denote the host plant from which the individual wasps originated. In the column ""HP", "Qv-F" indicates that the individual wasps are from Q. virginiana plants in sympatry. "Qv-T" indicates that the individual wasps are from Q. virginiana plants in allopatry. "Qg" designates that the individual wasps are from Q. geminata plants.

- The column labeled "Source..M." refers to the collection sites from which male wasps originated. The column labeled "Source..F." refers to the collection sites from which female wasps originated.

- Columns labeled as "C", "W", "Mo", "Mt", "Mc" represent the occurance of mating behavior: "contact","Wing buzz","male mount","copulation", "male courtship inlcuding wing buzz and mount". 0 means no such behavior occurred during 30 mininutes observation. 1 means such behavior occured at least once  during 30 mininutes observation.

- Columns labeled as "W.Latency", "Mo.Latency", "Mt.Latency" represent the proportion time it takes for the first observation of mating behavior: "Wing buzz","male mount","copulation".

- Columns labeled as "T.W", "T.Mo", "T.Mt", "T.Mc" represent the number of occurance of mating behavior: "contact","Wing buzz","male mount","copulation", "male courtship inlcuding wing buzz and mount".


