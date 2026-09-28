# FedGEE

**FedGEE** fits generalized estimating equations (GEE) across sites that
never pool patient data, such as hospitals in a health network. Sites share
only small $p \times p$ summary matrices. The package gives valid inference
even when there are only a few sites.

## Installation

```r
# install.packages("remotes")
remotes::install_github("soumikp/2026_fedGee")
```

## Features

- **Centralized** (`fedgee()`): a server sums the site summaries. The
  estimate equals pooled GEE.
- **Decentralized** (`decentralized_fedgee()`): no server. Sites average with
  neighbours over a gossip network (hub, ring, VISN-style or complete; see
  `build_weight_matrix()`).
- **Small-sample corrections in score space**: Kauermann–Carroll (`KC`,
  default), Mancl–DeRouen (`MD`) and Fay–Graubard (`FG`). All are computed
  from the site summaries alone.
- **Bell–McCaffrey degrees of freedom** (`df = "bm"`, default). They are
  computed from the site breads and adapt to unequal site sizes. KC with
  Bell–McCaffrey df equals CR2 with Satterthwaite df (clubSandwich) for
  linear models.
- **One fit, every variant**: `summary(fit, correction = "MD", df = "K-1")`
  switches variants without refitting.

## Example

```r
library(FedGEE)

data(ChickWeight)
cw <- as.data.frame(ChickWeight)
cw$site <- as.integer(cw$Chick) %% 12   # 12 mock sites
data_list <- split(cw, cw$site)

# Centralized
fit <- fedgee(data_list, weight ~ Time + Diet,
              family_obj = gaussian(), id_col = "Chick", verbose = FALSE)
fit                                      # KC + Bell-McCaffrey df
summary(fit, correction = "MD", df = "K-1")
confint(fit)

# Decentralized over a ring network
dfit <- decentralized_fedgee(data_list, weight ~ Time + Diet,
                             family_obj = gaussian(), id_col = "Chick",
                             structure = "ring", sandwich_level = "site",
                             correction = "KC",
                             L_beta = 150, L_S = 150, L_B = 150,
                             tol = 1e-6, verbose = FALSE)
dfit
```

## Notes

- The site-level sandwich has rank $\min(p, K - 1)$. With $K$ sites, keep
  $p$ well below $K$. `fedgee()` warns when the sandwich is singular.
- Each site estimates its own working correlation. The estimate equals
  pooled GEE exactly when every site uses the same correlation, which always
  holds for `corstr = "independence"`.
- Decentralized: disagreement between sites shrinks like $\rho^L$
  (`build_weight_matrix(...)$rho`). Too few rounds hurt the standard errors
  before they hurt the estimate.
