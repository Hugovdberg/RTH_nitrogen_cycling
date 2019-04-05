library(ReacTran)
library(marelac)
library(units)

aq_data <- read.table(
    'data/p4-data.dat',
    sep = "\t",
    col.names = c('distance', 'oxygen', 'ammonia')
)
plot(oxygen ~ distance, data = aq_data, type = 'l', lty = 1)
lines(ammonia ~ distance, data = aq_data, lty = 2)
