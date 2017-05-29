library(arulesSequences)
library(dplyr)

#loading data from file
diab.df <- read.csv("diab_trans.data", header=TRUE, stringsAsFactors = FALSE)

summary(diab.df)
head(diab.df)

#extract numeric codes from strings
diab.df$code <- as.integer(substr(diab.df$code, 4, 8))

#map meanings to codes
codes <- c()
codes["Regular insulin dose"]                     <- 33
codes["NPH insulin dose"]                         <- 34
codes["UltraLente insulin dose"]                  <- 35
#codes["Unspecified blood glucose measurement"]   <- 48 # same thing as 57
codes["Unspecified blood glucose measurement"]    <- 57
codes["Pre-breakfast blood glucose measurement"]  <- 58
codes["Post-breakfast blood glucose measurement"] <- 59
codes["Pre-lunch blood glucose measurement"]      <- 60
codes["Post-lunch blood glucose measurement"]     <- 61
codes["Pre-supper blood glucose measurement"]     <- 62
codes["Post-supper blood glucose measurement"]    <- 63
codes["Pre-snack blood glucose measurement"]      <- 64
codes["Hypoglycemic symptoms"]                    <- 65 # <------- IMPORTANT! we have to find rules with this result
codes["Typical meal ingestion"]                   <- 66
codes["More-than-usual meal ingestion"]           <- 67
codes["Less-than-usual meal ingestion"]           <- 68
codes["Typical exercise activity"]                <- 69
codes["More-than-usual exercise activity"]        <- 70
codes["Less-than-usual exercise activity"]        <- 71
codes["Unspecified special event"]                <- 72

#create reverse mapping as well
codeMeanings <- names(codes)
names(codeMeanings) <- codes

# unify two values meaning the same meaning ("Unspecified blood glucose measurement")
diab.df$code[diab.df$code==48] <- 57

#filter unknown codes
diab.df <- diab.df[ diab.df$code %in% codes, ]

# add code meaning
diab.df$code_meaning = factor(diab.df$code, levels = codes, labels = codeMeanings)
head(diab.df$code_meaning)
diab.df$code_meaning <- as.character(diab.df$code_meaning)

#Some kinds of events we have to deal separately because value field means something different

##########    glucose measurment related records    ##########
diab.glucose <- diab.df[ (unlist(diab.df$code >= 57) & unlist(diab.df$code <= 64)),]
#filter deadly low values - (possibly an error in data)
diab.glucose <- diab.glucose[diab.glucose$value>10,]
#split values into 5 ranges - very low, low, medium, high, very high
diab.glucose$value <- cut(diab.glucose$value, breaks = 10)

#create event description for glucose measurments
diab.glucose$eventDescr <- paste0( as.character(diab.glucose$code_meaning),' value in range ', diab.glucose$value)
head(diab.glucose)

##########     insulin doses related records        ##########
diab.dose <- diab.df[ (unlist(diab.df$code >= 33) & unlist(diab.df$code <= 35)),]
#filter out deadly high doses of insulin or negative values - (possibly an error in data)
diab.dose <-(diab.dose[ (unlist(diab.dose$value <= 80) & unlist(diab.dose$value > 0)  ),])
#split values into 5 ranges - very low, low, medium, high, very high
diab.dose$value <- cut(diab.dose$value, breaks = 5)
#create event description for insulin doses
diab.dose$eventDescr <- paste0( as.character(diab.dose$code_meaning),', units count in range ', diab.dose$value)

head(diab.dose)
summary(diab.dose)

##########     events without meaningfull value     ##########
diab.other <- diab.df[ (unlist(diab.df$code >= 65) & unlist(diab.df$code <= 72)),]
diab.other$eventDescr <- paste0( as.character(diab.other$code_meaning) )

head(diab.other)
summary(diab.other)

##########   Fuse data together to single table    ##########
diab.df <- rbind(diab.glucose, diab.dose, diab.other)
#filter not-NA values
diab.df<- diab.df[unlist(diab.df$eventDescr != "NA value in range NA"),]
diab.df<- diab.df[unlist(!is.na(diab.df$patient_id)),]
diab.df<- diab.df[unlist(!is.na(diab.df$time_sek)),]

unique(diab.df$eventDescr)
head(diab.df)
summary(diab.df)

##########     Transform data to sequence form     ##########

#drop no longer needed collumns.
diab.df$code <- NULL
diab.df$code_meaning <- NULL
diab.df$value <- NULL

head(diab.df)
summary(diab.df)

#sort by sequence and event ID (Needed by spade impl)
diab.df <- diab.df[with(diab.df, order(patient_id,time_sek)), ]

head(diab.df)

#write data to CSV
write.table(diab.df, "diab_trans_processed.data", sep = "@", quote = FALSE, row.names = FALSE, col.names = FALSE )

#read as transaction sequence form
diabSeq <- read_baskets(con = "diab_trans_processed.data", sep ="@", info = c("sequenceID","eventID"))
summary (diabSeq )
head (diabSeq )

##################    Detect sequences    ###################

#params of sequence detection
seqParam <- new ("SPparameter",support = 0.03, maxsize = 5, mingap=1, maxgap =14400, maxlen = 8)

#detect common sequence patterns using SPADE algorithm
patSeq <- cspade(diabSeq, seqParam, control = list(verbose = TRUE, tidLists = TRUE, summary= TRUE))
inspect(patSeq)

#rule generation basing on common sequences
seqRules = ruleInduction(patSeq,confidence = 0.8)

#all sequences
allSeq <- c(rhs(seqRules),lhs(seqRules))
allSeq <- unique(allSeq)
inspect(allSeq)
str(allSeq)

inspect(lhs(seqRules))
inspect(rhs(seqRules))

##################    Select sequences    ###################
#select interesting sequences - the one with hypoglycemic symptoms as a results
rulesI = subset(seqRules, rhs(seqRules) %in% c("Hypoglycemic symptoms"))
inspect(rulesI)
View(as(rulesI,"data.frame"))
