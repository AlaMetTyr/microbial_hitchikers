library(glmmTMB)
library(ggplot2)

data <- read.delim("clipboard")

##modelling state versus count
model_state <- glmmTMB(
  Count ~ State + (1 | Reps),
  offset = log(Dilution),
  family = nbinom2,
  data = data
)
summary(model_state)


##modelling method with state (clean)
clean_only <- subset(data, State == "clean")

model_method <- glmmTMB(
  Count ~ Method + (1 | Reps),
  offset = log(Dilution),
  family = nbinom2,
  data = clean_only
)
summary(model_method)


#modelling site versus state (is one worse than others)
dirty_only <- subset(data, State == "dirty")

model_site <- glmmTMB(
  Count ~ Site + (1 | Reps),
  offset = log(Dilution),
  family = nbinom2,
  data = dirty_only
)
summary(model_site)

##maybe make a graph
library(ggplot2)
ggplot(data, aes(x = Method, y = Count)) +
  geom_boxplot() +
  labs(title = "Counts by State")
