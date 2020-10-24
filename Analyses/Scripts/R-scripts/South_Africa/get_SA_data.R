# get South African data to use with scripts
# treedater_gag.R; treedater_pol.R and treedater_env.R
library(lubridate)

pangea_data <- read.csv("~/Box Sync/pangea2018/data/PANGEA_Extract_Volz_2019-04-16.csv",
                        na.strings=c("","NA"))



# get data from South Africa sequences
sa <- subset(pangea_data, main_cohort_id == "South Africa" | geo_country == "SouthAfrica")
sa <- subset(sa, gender != "Unknown")
sa$pangea_id <- as.character(sa$pangea_id)
sa$pangea_id <- as.factor(sa$pangea_id)

sa$geo_country <- as.character(sa$geo_country)
sa$geo_country <- as.factor(sa$geo_country)

sa$gender <- as.character(sa$gender)
sa$gender <- as.factor(sa$gender)

# get pangea data that was used as CGR sequences
load( system.file( "data/SA/SA_pangea_cgr.rda", package = 'pangea'))


# get sample_date for each South African sequence
sa_dates <- sa[c("pangea_id", "geo_country", "sample_date", "shiver_consensus")]

# check if there is any NA in the sample dates
sum(is.na(sa_dates$sample_date)) # no NA in sample dates

# convert dates to year-month-day format
sa_dates["sa_ymd"] <- ymd(sa_dates$sample_date)

# convert dates to decimal format (because treedater used dates in decimal
# format)
sa_dates["decimal"] <- decimal_date(sa_dates$sa_ymd)
sa_dates$decimal <- as.factor(sa_dates$decimal)
