# Code for _How good is 85%? A survey tool to connect classifier evaluation to acceptability of accuracy_

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img
alt="This work is licensed under a Creative Commons Attribution 4.0 International License" 
title="This work is licensed under a Creative Commons Attribution 4.0 International License" 
src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a> _Matthew Kay ([mjskay@uw.edu]
(mailto:mjskay@uw.edu)), Shwetak N. Patel ([shwetak@uw.edu]
(mailto:shwetak@uw.edu)), and Julie A. Kientz ([jkientz@uw.edu]
(mailto:jkientz@uw.edu))_

This repository contains analysis code from: 

Kay, Matthew, Patel, Shwetak N., and Kientz, Julie A. How Good is 85%? A Survey 
Tool to Connect Classifier Evaluation to Acceptability of Accuracy. _CHI 2015_
(upcoming). http://dx.doi.org/10.1145/2702123.2702603

It is intended to allow others to adopt our tool for generating surveys and
modelling acceptability of accuracy. It is currently a work-in-progress. If you
have any questions, please email Matthew Kay (above). Also, if you've done 
something cool with this work, we'd love to hear from you!

## Libraries needed for these analyses




```r
library(runjags)
library(tidybayes)
library(coda)
library(lme4)
library(plyr)
library(dplyr)
library(ggplot2)
library(pander)
library(GGally)

memory.limit(8000)
```

```
## [1] 8000
```

Plus some helper functions in [util.R](util.R):


```r
source("util.R")
```

## Example model

### Fitting the model
Example code for fitting a model can be found in [src/application_ui-regression.R](src/application_ui-regression.R). Since fitting 
the model takes some time, for the purposes of this example we will just load the fitted model, 
which has been saved in `src/output/acceptability_ui-model-small-final.RData`:


```r
load(file="output/acceptability_ui-model-stan-final.RData")
fit = apply_prototypes(fit, df)
```

### Parameter estimates

Let's plot some posterior parameter estimates. First, we extract the 
sample estimates for `b0`, `b`, `p`, and `alpha`:


```r
params = extract_samples(fit, cbind(b, alpha_mu, alpha_v, p)[application])
thresholds = extract_samples(fit, b0[application,acceptable])
sigmas = extract_samples(fit, cbind(sigma_b, sigma_b0)[])
```

This gives us a pretty simple table of estimates. We can look at the
first couple of entries to get an idea of its structure:


```r
head(params)
```


|  .sample  |    application     |  alpha_mu  |  alpha_v  |   b   |   p    |
|:---------:|:------------------:|:----------:|:---------:|:-----:|:------:|
|     1     |    alarm_police    |   0.4757   |   3.678   | 17.82 | 0.3606 |
|     1     | alarm_text_message |   0.3576   |   1.944   | 18.55 | -2.301 |
|     1     |    electricity     |   0.6023   |   4.433   | 15.82 | -2.033 |
|     1     |      location      |   0.5265   |   1.919   | 15.69 | -2.733 |
|     2     |    alarm_police    |   0.5009   |   2.25    | 16.95 | -0.286 |
|     2     | alarm_text_message |   0.4082   |   1.186   | 17.59 | -1.747 |

Now we'll plot each parameter in turn, with some useful reference lines:


```r
ggposterior(params, aes(x=application, y=b)) +
    geom_hline(yintercept=0, lty="dashed")
```

![plot of chunk params_plot](figure/params_plot-1.svg) 

```r
ggposterior(thresholds, aes(x=acceptable, y=b0)) +
    geom_hline(yintercept=0, lty="dashed") +
    facet_wrap(~application) 
```

![plot of chunk params_plot](figure/params_plot-2.svg) 

```r
ggposterior(params, aes(x=application, y=p)) +
    geom_hline(yintercept=-1, lty="dashed")
```

![plot of chunk params_plot](figure/params_plot-3.svg) 

```r
ggposterior(params, aes(x=application, y=alpha_mu)) +
    geom_hline(yintercept=0.5, lty="dashed") +
    ylim(0, 1)
```

![plot of chunk params_plot](figure/params_plot-4.svg) 

```r
ggposterior(params, aes(x=application, y=alpha_v))
```

![plot of chunk params_plot](figure/params_plot-5.svg) 

```r
params %<>% mutate(alpha_var = alpha_mu * (1 - alpha_mu)/(alpha_v + 1))
```

```
## Error in eval(expr, envir, enclos): could not find function "%<>%"
```

```r
ggposterior(params, aes(x=application, y=alpha_var))
```

```
## Error in eval(expr, envir, enclos): object 'alpha_var' not found
```

```r
ggposterior(sigmas, aes(x=1, y=sigma_b))
```

![plot of chunk params_plot](figure/params_plot-6.svg) 

```r
ggposterior(sigmas, aes(x=1, y=sigma_b0))
```

![plot of chunk params_plot](figure/params_plot-7.svg) 

```r
ggpairs(log(cbind(select(sigmas, -.sample))),
    lower = list(continuous = "density"))
```

![plot of chunk params_plot](figure/params_plot-8.svg) 

### Differences in alpha

We can also examine the posterior difference in alpha between each condition. For 
every pair of applications, let's get the posterior distribution of their 
difference in alpha:


```r
alpha_comparisons = compare_levels(params, alpha_mu, by=application)
```

Which looks like this:


```r
head(alpha_comparisons)
```


|  .sample  |            application            |  alpha_mu  |
|:---------:|:---------------------------------:|:----------:|
|     1     | alarm_text_message - alarm_police |  -0.1181   |
|     2     | alarm_text_message - alarm_police |  -0.09262  |
|     3     | alarm_text_message - alarm_police |  -0.2586   |
|     4     | alarm_text_message - alarm_police |  -0.2536   |
|     5     | alarm_text_message - alarm_police |  -0.05727  |
|     6     | alarm_text_message - alarm_police |  -0.02978  |

Finally, we can plot the estimated differences:


```r
ggposterior(alpha_comparisons, aes(x=application, y=alpha_mu)) + 
    geom_hline(yintercept=0, lty="dashed") +
    ylim(-0.6, 0.6)
```

![plot of chunk alpha_comparison_plot](figure/alpha_comparison_plot-1.svg) 

## Citing this work

Please cite the CHI paper above.

## Problems

Should you encounter any issues with this code, contact Matthew Kay
(<mjskay@uw.edu>). If you have found a bug, please file it [here](https://github.com/mjskay/acceptability-of-accuracy/issues/new) with minimal code to reproduce
the issue.


