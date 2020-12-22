####################
# Project overview #
####################

The following files are associated with the project "A global meta-analysis of temperature effects on 
marine fishes’ digestion across trophic groups" and are intended to allow any interested parties to 
access the data and code used to produce the analyses in that paper.

All data files are accessible through Dryad.
All data files, analysis files and .rds files are accessible through Github at https://github.com/nicoleknight0/meta-analysis-fish-digestion

Contact nicole.knight@mail.mcgill.ca with any questions.


##################
# Analysis files #
##################

### data_analysis.R
-R script that provides code that can be used to reproduce all final models, 
analyses, figures and tables in the published main text and Appendix 2.  

### sessionInfo.txt
-Text file that provides R version and package versions used to conduct the final models and analyses

##############
# .rds files # 
##############

All .rds files are records of the models generated from the data_analysis.R script using set.seed(438).
Parameter estimates used in the figures and tables included in the main text and Appendix 2 were 
extracted from models built with that specific set.seed() and run in the same order as presented in the
script.  Running the models without set.seed(438) may slightly change parameter estimates. 

##############
# Data files #
##############

### data_absorption_efficiency.csv

-This file contains all estimates of marine fishes' absorption efficiency and ancillary data taken from the literature
that were used in our analyses.

### data_gut_length.csv

-This file contains all estimates of marine fishes' gut length and ancillary data taken from the literature that 
were used in our analyses.

### data_gut_passage.csv

-This file contains all estimates of marine fishes' gut passage time and ancillary data taken from the literature that 
were used in our analyses.

### data_nutrients.csv

-This file contains all estimates of the nutrient concentrations of food fed to fishes in studies 
included in our other analyses and ancillary data. 

### model_weights.csv

-This file contains a record of all model weights used to calculate the average and standard deviations
of the model weights reported in the main text and Appendix 2. 

############################
# Data file column headers #
############################

#Across all files

species	- > study species 				
study.id -> arbitrary numeric ID for study
dataset.id -> arbitrary numeric ID for dataset (separate experiments within study)
wild.lab -> data collected from individuals in wild or lab?
order -> phylogenetic order of study species
family -> phylogenetic family of study species
genus -> phylogenetic genus of study species
trophic.group -> general classification: carnivore, herbivore or omnivore
latitude -> approximate latitude that the fish were collected from
longitude -> approximate longitude that the fish were collected from
study.location -> location that the study was conducted 
diet -> more specific classification: type of plant or animal material consumed
diet.details -> further details on diet
temperature.c -> temperature of collection site or experimental setup.  When not reported by the original study, collected from NOAA (see methods)
temperature.variation -> reported variability in temperature at collection site or in experimental setup
age.class -> juvenile, subadult, adult, etc.
bodylength.metric -> standard, total or fork length (entire body)
ln.bodylength -> natural logarithm of average body length
bodylength.avg.mm -> average body length in millimeters (or midpoint of length range, these instances are noted in the notes column)
bodylength.min.mm -> minimum body length in millimeters
bodylength.max.mm -> maximum body length in millimeters
bodylength.sd.mm -> standard deviation of body length in millimeters
bodylength.se.mm -> standard error of body length in millimeters
bodylength.n -> sample size for body length
mass.avg.g -> average body mass in grams
mass.sd.g -> standard deviation of body mass in grams
mass.min.g -> minimum body mass in grams
mass.max.g -> maximum body mass in grams 
mass.n -> sample size for body mass
reference -> full reference for study that provided the original data 
notes -> additional details on the dataset or data collection
more notes -> additional details on the dataset or data collection

#Specific to data_absorption_efficiency.csv

method -> for absorption efficiency, marker method or total collection method (see Appendix S2 of publication)
absorption.efficiency.n -> sample size for absorption efficiency
component -> component absorbed (e.g., energy or protein or organic)
absorption.efficiency -> proportion of component absorbed (0 = none of the component absorbed, 1 = all of the component absorbed)

#Specific to data_gut_length.csv

preferred -> for species with multiple gut length reports, poorer-quality reports (e.g., did not report temperature) or reports of intestine-only length were excluded (coded N) and the script randomly selects from the remaining reports (coded Y)
type -> did the authors measure only the intestine or the entire gut
diet.reference -> if the species' diet classification was not reported in the original study, the study from which we collected its diet classification
sample.period -> the period over which fish were collected for gut length analysis
gutlength.mean.mm -> average gut length in millimeters
ln.gutlength.mean -> natural logarithm of gut length
gutlength.sd -> standard deviation of gut length
relativelength.min -> minimum reported (gut length) / (body length)
relativelength.max -> maximum reported (gut length) / (body length)
relativelength.mean -> mean reported reported (gut length) / (body length)
relativelength.sd -> standard deviation of reported (gut length) / (body length)
relativelength.n -> sample size of reported (gut length) / (body length)
ln.relativelength -> natural logarithm of reported (gut length) / (body length)

### data_gut_passage.csv

diet2 <- diet classification used to produce figure 2.  Distinguishes known fermenters and non-fermenters in the dataset
passage.time.h <- time required after feeding for a parcel of food to pass through the entire digestive tract in hours
ln.passage.time	<- natural logarithm of gut passage time
passage.time.sd	<- standard deviation of gut passage time
passage.time.n <- gut passage time sample size

### data_nutrients.csv
diet.type <- food item plant or animal
protein <- protein in food as a proportion of dry mass
lipid <- lipid in food as a proportion of dry mass
carbohydrate <- carbohydrates in food as a proportion of dry mass
nitrogen <- nitrogen in food as a proportion of dry mass
carbon <- carbon in food as a proportion of dry mass
organic <- organic in food as a proportion of dry mass
ash <- ash in food as a proportion of dry mass
energy.reported <- energy per unit dry mass as originally reported
energy.reported.units <- units used in original energy report
energy.kJg <- kJ energy per gram of dry mass (converted from original report to standardize units)
