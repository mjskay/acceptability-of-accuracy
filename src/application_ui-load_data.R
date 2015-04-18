## Loads data for the application and UI face validity experiment
# 
# Author: mjskay
###############################################################################

df = read.csv("data/application_and_ui.csv")

#code likerts as ordinal factors
for (column in c("acceptable", "useful", "would_use")) {
	df[[column]] = ordered(df[[column]], 
        levels=c("Extremely unlikely", "Quite unlikely", "Slightly unlikely", "Neither", "Slightly likely", "Quite likely", "Extremely likely"))
}
