require(data.table)
require(reshape)
require(ggplot2)

#Run install.packages("data.table", "reshape", "ggplot2")
#if you are missing the above packages.

res <- data.table(read.csv("../io/out.csv"))
mres <- data.table(melt(res, id.vars='iteration'))
qplot(iteration, value, data = mres[variable != "likelihood"], geom='step') +
    facet_grid(~variable, ncol=2, scales="free_y", labeller = label_parsed) +
  opts(title=expression("Heston parameters for IBM as calculated by EKF"))

