### Helper functions for making nice plots

#function for use with stat_summary that returns the median and
#highest density interval of the data
median_hdi = function(x, ...) {
    coda:::HPDinterval.mcmc(x, ...) %>% 
        data.frame() %>% 
        select(ymin = lower, ymax = upper) %>% 
        cbind(y = median(x, ...))
}

#simple posterior violin plot
ggposterior = function(.data, .aes) {
    ggplot(
            .data,
            .aes
        ) + 
        geom_violin(linetype=0, fill="skyblue") + 
        stat_summary(fun.data=median_hdi) +
        coord_flip() 
}
