install.packages("corrplot")
install.packages("ggplot2")
install.packages("ggpubr")

library(corrplot)
library(ggplot2)
library(ggpubr)

# Load data
td <- read.csv( "data of PPI network analysis.csv")

# correlation analysis
tdc <- cor (td, method="pearson")

# Data visualization

corrplot(tdc, method = "circle", type = "upper",
         tl.col = "black", tl.cex = 0.7, tl.srt = 60)

corrplot(tdc, method = "number", type = "lower", 
         tl.col = "n", tl.cex = 0.6, tl.pos = "n",order = 'AOE',
         add = T)

# Save data
write.csv(tdc,"correlation analysis.csv")





