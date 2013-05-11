require(data.table)
require(reshape)
require(ggplot2)

res <- data.table(read.csv("../io/out.csv"))
mres <- data.table(melt(res, id.vars='iteration'))
qplot(iteration, value, data = mres[variable != "likelihood"], geom='step') +
    facet_wrap(~variable, ncol=2, scales="free_y")
