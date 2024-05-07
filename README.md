# Multiple outcome Synthetic Difference-in-Differences

Allows to use multiple outcomes in Synthetic Diff-in-Diffs estimation. 
Implements the synthetic difference-in-differences estimation proposed by [Chava and Reguly (2024): The Long-Run Stock Market Performance of Mergers and Acquisitions](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4742401)
The package extends [Arkhangelsky et al (2021)](https://www.aeaweb.org/articles?id=10.1257/aer.20190159) to multiple outcomes, while optimizing the weights over multiple variables.

## How is it different from classical Synthetic Diff-in-Diffs

Arkhangelsky et al (2021) introduces synthetic control that allows time and unit fixed effects while creating ''synthetic control'' for the (expected) treated unit.

Other Synthetic Control methods allows to use ''covariates'' that they use to calculate the weights. 
Note however these variables are unaffected by the treatment. Nice paper on the issue is by [Botosaru and Ferman (2019)](https://academic.oup.com/ectj/article-abstract/22/2/117/5303851?login=false).

Here we allow to estimate the weights for the control donor pool by using multiple outcome that are affected by the treatment. We also allow pre-event time specific weights that reduces the bias in the estimator.

$$ \hat Y_{it}(1) = \sum_{i'} \omega_{i'} Y_{i't} - \hat \delta _{pre} $$

$$ \hat \delta_{pre} = \sum_{t=1}^{T_{pre}} \left( \hat\lambda_{t}Y_{i,t} - \sum_{i'} \hat \omega_{i'} \hat\lambda_{t} Y_{i't}\right) - \hat\lambda_{0,i} $$

where the individual unit-weights $\omega_{i'}$ and the pre-treatment time weigths $\lambda_{t}$ are estimated by using multiple outcomes. 
Their optimal values are the solution of the following optimization problem:

$$\hat\omega_{i'} = \left( \hat\omega_0^K, \hat\omega_{i'} \right) = arg\min_{\omega_0\in \mathbb{R}, \omega \in \Omega} \sum_{t=1}^{T_{pre}} 
\left[ \Big(\omega_0 + \sum_{i'} \omega_{i'} Y_{i',t} - Y_{i,t} \Big)^{'} \hat \Sigma_t^{-1} \Big( \omega_0 + \sum_{i'} \omega_{i'} Y_{i',t} - Y_{i,t} \Big) \right] + \zeta^2 T_{pre} \sum_{i'}\omega^2_{i'}$$

and 

$$ \left( \hat\lambda_0, \hat\lambda \right) = arg\min_{\lambda_0 \in \mathbb{R}, \lambda \in \Lambda}
    \sum_{i'} \Bigg[ \left(\lambda_0+\sum_{t=1}^{T_{pre}}\lambda_{t} Y_{i',t} - \frac{1}{T_{post}}\sum_{t=0}^{T_{post}} Y_{i',t}\right)^{'}
    \hat{\Sigma}^{-1} \left(\lambda_0+\sum_{t=1}^{T_{pre}}\lambda_{t} Y_{i',t} - \frac{1}{T_{post}}\sum_{t=0}^{T_{post}} Y_{i',t}\right) \Bigg]$$

where $Y$ is $J \times T$, $i$ is the treated units and $i'$ are the controls and during the optimization we constraint both weights to sums up to 1. See more in [Chava and Reguly (2024)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4742401)

Intuitively it uses information from $J$ outcome and calculate unit specific weigths $\omega_{-i}$ that minimize the Mahalanobis distance 
from the treated (average) unit, where the scaler is coming from time specific var-covar matrices. 
Similarly to Arkhangelsky et al (2021) we add a (ridge) penalty term that shifts weights toward zero.

For the time specific weights we use the same approach but scale with unit-specific var-covar matrices.

# Usage

## Installation

To install this package, first ensure that devtools is installed with

```
install.packages("devtools")
```

then install the package from GitHub with

```
devtools::install_github("regulyagoston/multisynthdid")
```

## Basic usage

To get started, use a panel dataset with an outcome measure, a treatment indicator, a unit indicator, a time variable, and an intervention time t_int. Then run

```
res <- multi_sdid( outcome_1 + outcome_2 ~ 1, 'treatment', 'unit_id', 'time_id', data )
```

## Example with German Reunification Data

A simple example is to estimate the effect of German Reunification.

```
library(multisynthdid)
```

Get the imputed data for german reunification.

```
data("german_reunification")
```

One can easily run the estimation via the main function:

```
res <- multi_sdid( gdp + infrate + trade ~ 1, 'W', 'country', 'year', german_reunification )
```

The summary command will summarise the main estimates

```
summary.multi_sdid(res)
```

To check the outcome variables and its synthetic counterparts visual investigation has been the norm with synthetic controls.

```
plot_outcomes( res )
```

Finally one can check the treatment/gap plots with ot without placebo effects for inference.

```
plot_gaps( res, placebo = T )
```


