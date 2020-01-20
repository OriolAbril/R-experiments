# Use mombf to select important variables and to predict silica concentration
# on Kaggle dataset "Quality Prediction in a Mining Process"
library("mombf")
data <- read.csv("MiningProcess_Flotation_Plant_Database.csv")
data <- data.frame(data[, !names(data) %in% c("date","X..Iron.Concentrate")])
data[] <- lapply(
  data,
  function(x) {
    if(is.factor(x)) as.numeric(as.character(gsub(",", ".", x, fixed=TRUE))) else x
  }
)

p <- length(names(data))
formula_text <- paste("X..Silica.Concentrate", paste0(names(data)[1:p-1], collapse=" + "), sep=" ~ ")

y <- data[,p]
X <- as.matrix(data)
X[,p] <- 1

# fit <- modelSelection(as.formula(formula_text), data=data, )
fit <- modelSelection(y=y, x=X, niter=20000)
postProb(fit)


