
library(glmmTMB)
library(ggplot2)

data <- read.delim("clipboard")

# Log-transform dilution
data$LogDilution <- log(data$Dilution)

##model nested
model <- glmmTMB(
  Count ~ State + Method + Site + State:Method + State:Site + offset(LogDilution) + (1 | Reps),
  family = nbinom2,
  data = data
)
##or?
# Model with nested structure, including LogDilution as offset and random intercept for Reps
model_nested <- glmmTMB(
  Count ~ State * Method + State * Site + offset(LogDilution) + (1 | Reps),
  family = nbinom2,
  data = data
)
summary(model_nested)

## model simple
model <- glmmTMB(Count ~ State + offset(LogDilution) + (1 | Reps), 
                 family = nbinom2, data = data)

##model including log10 as a random effect?
model <- glmmTMB(Count ~ State + LogDilution + (1 | Reps), 
                 family = nbinom2, data = data)
summary(model)


#model summary
summary(model)

table(data$State)
boxplot(Count ~ State, data = data)

ggplot(data, aes(x = State, y = Count)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution of Counts by State")

ggplot(data, aes(x = Count, fill = State)) +
  geom_histogram(binwidth = 1, position = "dodge") +
  theme_minimal() +
  labs(title = "Histogram of Counts by State")

boxplot.stats(data$Count[data$State == "dirty"])

# Plot residuals
plot(resid(model))


###now llets try fix the model
##############################




model_nested <- glmmTMB(
  Count ~ State * Method + (1 | State:Site) + offset(LogDilution) + (1 | Reps),
  family = nbinom2,
  data = data
)

#dilution as offset
model_offset <- glmmTMB(
  Count ~ State + (1 | Reps),
  offset = LogDilution,
  family = nbinom2,
  data = data
)

#dilutiona s fixed effect
model_fixed <- glmmTMB(
  Count ~ State + LogDilution + (1 | Reps),
  family = nbinom2,
  data = data
)

##compare logdiultion as offset or fixed effect
AIC(model_offset, model_fixed)

##“We used dilution as an offset, since it directly scales the concentration of organisms and thus the expected count.”
##If reviewers or collaborators ask, you can also say:
##“We compared both offset and fixed-effect approaches. While the fixed-effect model had a slightly lower AIC, we chose the offset because the effect of dilution is known and expected to scale linearly on the log scale.”
model <- glmmTMB(
  Count ~ State + (1 | Reps),
  offset = log(Dilution),
  family = nbinom2,
  data = data
)
summary(model)

##model including testing cleaning methods
model_method <- glmmTMB(
  Count ~ State + Method + (1 | Reps),
  offset = log(Dilution),
  family = nbinom2,
  data = data
)
summary(model_method)

