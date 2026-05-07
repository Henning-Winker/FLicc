# FLicc

**FLicc: Length-interval catch-curve per-recruit analysis in FLR**

`FLicc` is an user-friendly TMB implementation in the `FLR` and a merged generalization of both the multi-gear length-interval catch-curve framework developed in [`fishblicc`](https://github.com/PaulAHMedley/fishblicc) ([Medley 2025](https://academic.oup.com/icesjms/article/82/12/fsaf224/8404907?guestAccessKey=)) and the length-based spawning potential ratio estimation in [`LBSRP`](https://github.com/adrianhordyk/LBSPR) (Hordyk et al. [2015](https://academic.oup.com/icesjms/article/72/1/217/815772?guestAccessKey=.). [2016](https://cdnsciencepub.com/doi/full/10.1139/cjfas-2015-0422)).

`FLicc` retains the biological core of either of the original length-based equilibrium models, but adds a fast penalized-likelihood workflow in Template Model Builder (TMB), multi-year fitting, FLR integration, and a growing set of equilibrium and indicator tools. The original multi-gear formulation was designed to estimate mortality-at-length, selectivity, spawning potential ratio (SPR), yield-per-recruit from length compositions grouped by gear, and parallel computing utilities. 

An introductory User Guide to Multi-Gear Length-Interval Catch-Curve Analysis in FLR for the `FLicc` package is available [here](https://henning-winker.github.io/FLicc/FLicc_intro_vignette.html). 

An overview with illustrations of the FLicc utilities is provided [here](https://github.com/Henning-Winker/FLicc/blob/main/docs/FLicc_Intro.pdf).

---

## Motivation

The original motivation is simple: many real fisheries are **multi-gear**, but many length-based methods, such as `LBSRP` assume a single gear with logistic selectivity.

Many length-based methods are built around a single gear and logistic selectivity. That can work well in simple fisheries, but it is often a poor description of real fleets, especially when gears differ strongly in size selectivity or when one or more gears are dome-shaped. Medley’s multi-gear framework was developed specifically to address this problem by jointly fitting gear-specific length compositions and relative catches from multiple fleets. 

This matters because different gears contain different information:

- **Fishing impact** is informed by the relative catch taken by each gear.
- **Population size structure** is informed by the combined shape of the observed length distributions across gears.
- **Selectivity shape** is better identified when gears overlap in the stock but differ in how they sample length classes.

In practice, combining gears can help separate mortality from selectivity, improve estimation of stock structure, and reduce the risk that dome-shaped selectivity is misread as a pure mortality signal. Medley’s simulation results showed that multiple gears with different selectivities can provide better SPR estimates than relying on a single gear, especially when at least one gear still samples larger fish. 

---

## Relationship to LBSPR and per-recruit models

`FLicc` is closely related to the **length-based spawning potential ratio (LBSPR)** framework developed by Hordyk et al. (2015). Both approaches are based on a **length-structured per-recruit model**, where population structure and SPR are determined primarily by the ratio:

$$
\frac{Z}{K} = \frac{M + F}{K}
$$

A key result of this framework is that abundance-at-length can be expressed recursively as:

$$
N_L \propto (L_\infty - L)^{Z/K}
$$

This allows SPR and fishing mortality to be estimated directly from length composition data without requiring age-structured models.

### Reproducibility

Under the special case of:

- a **single fleet**  
- **logistic selectivity**  
- a **multinomial likelihood**  

`FLicc` reproduces LBSPR-type inference. In this sense, `FLicc` provides a reproducible and extensible implementation of the LBSPR framework within an FLR/TMB workflow.

### Extensions in FLicc

Relative to LBSPR, `FLicc` extends the framework by:

- supporting **multiple gears with distinct selectivity patterns**,  
- allowing **flexible and dome-shaped selectivity**,  
- separating **population and observation models**,  
- enabling **multi-year estimation**, and  
- providing alternative likelihoods for overdispersed data.

These extensions improve identifiability of fishing mortality and selectivity, particularly in mixed-gear fisheries.



---

## Why multi-gear fitting helps

A central advantage of `FLicc` is that multi-gear data can inform different parts of the inference problem at the same time.

If only one gear is fit, the observed descending limb of the catch-at-length distribution can reflect both:
1. real population depletion with size, and
2. the descending side of selectivity.


Estimation of fishing mortality from a single catch-at-age or length composition can be weakly identified when selectivity is dome-shaped. In such cases, the descending limb of the observed distribution reflects a combination of total mortality and declining selectivity at older ages or larger sizes. This induces confounding between fishing mortality and selectivity, such that different combinations of these processes can produce similar observed patterns.

In a multi-gear setting, each fleet samples the population through a distinct selectivity pattern. When these patterns differ in their coverage of the size or age range—particularly when at least one gear provides information on larger or older individuals—the combined data help disentangle selectivity from mortality.

The key mechanisms are:

- **Contrast in selectivity shapes**: Differences in ascending and descending limbs across gears provide independent information on population structure.
- **Coverage of large individuals**: Gears that retain larger fish reduce ambiguity associated with the descending limb of dome-shaped selectivity.
- **Shared population dynamics**: All gears observe the same underlying population, imposing consistency across inferred mortality processes.

As a result, the joint likelihood across gears constrains the solution space, reducing the risk of spurious selectivity patterns and improving the robustness of fishing mortality estimates.

---

## Core model idea

The model assumes an approximate steady state and tracks survival through sequential length intervals under growth and mortality. In the underlying framework, survival to length interval $n$ is computed across intervals and integrated over variability in asymptotic length $L_\infty$, allowing growth variability to be carried into the expected length composition. Numerical integration is then used to obtain relative abundance-at-length. :contentReference[oaicite:4]{index=4}

Catch in length interval $i$ is given by:

$$
c_i = \frac{F_i}{Z_i}(S_i - S_{i+1})
$$

where $S_i$ is survival to the lower bound of interval $i$, $F_i$ is fishing mortality at length, and $Z_i = M_i + F_i$ is total mortality at length. For gear $j$, catch per recruit is:

$$
C_j = \sum_{i=1}^{N} \frac{F_{ij}}{Z_i}(S_i - S_{i+1})
$$

and the relative catch share of gear $j$ is:

$$
p_j = \frac{C_j}{\sum_{g=1}^{G} C_g}
$$

so the relative catches provide information on the relative fishing pressure exerted by each gear. 

Total mortality at length is modeled as the sum of natural mortality and fishing mortality across gears:

$$
Z_i = M_i + \sum_j F_j \sum_k w_k \, sel_{ijk}
$$

where $F_j$ is the apical fishing mortality for gear $j$, $sel_{ijk}$ is the selectivity component at length for gear $j$, and $w_k$ are mixture weights when multiple selectivity components are used. This is the key decomposition that allows `FLicc` to separate gear impact from information on stock size structure. :contentReference[oaicite:6]{index=6}

---



---

## Population and observation models

`FLicc` separates the model into:

- a **population model**, which generates the expected length composition, and  
- an **observation model**, which links observed data to model expectations.

### Default configuration

```r
settings = list(pop_model = "gtg", obs_model = "mn")
```

This configuration is typically **5–10× faster** than the original fishblicc formulation and is recommended for most applications.

---

### Population models

#### GTG (default)

The **growth-type-group (GTG)** formulation approximates growth variability and survival efficiently:

$$
\pi_i = \frac{N_i}{\sum_k N_k}
$$

with recursive survival:

$$
N_{i+1} = N_i \exp(-Z_i \Delta t_i)
$$

This formulation captures the key structure of LBSPR-style models while avoiding costly numerical integration.

---

#### Gamma (fishblicc-style)

The original fishblicc formulation integrates over variability in asymptotic length:

$$
N_i = \int f(L_\infty) \, S_i(L_\infty) \, dL_\infty
$$

This approach is more flexible but computationally slower due to numerical integration.

---

### Observation models

Let $y_i$ be observed counts and $\pi_i$ expected proportions.

#### Multinomial (default)

$$
\mathbf{y} \sim \text{Multinomial}(n, \boldsymbol{\pi})
$$

Fast, stable, and consistent with LBSPR.

---

#### Dirichlet–multinomial

$$
\mathbf{y} \sim \text{DM}(n, \boldsymbol{\pi}, \theta)
$$

Allows for overdispersion in length composition data.

---

#### Negative binomial

$$
y_i \sim \text{NB}(\mu_i, \phi), \quad \mu_i = n \pi_i
$$

Flexible variance structure but slower and less stable.

---

### Practical guidance

- **Default (recommended):** GTG + multinomial  
- **fishblicc replication:** gamma + negative binomial  
- **Overdispersed data:** Dirichlet–multinomial  



## Natural mortality at length

`FLicc` supports alternative formulations for natural mortality as a function of length, $M(l)$, allowing users to explore different assumptions about size-dependent mortality.

Let $M_{ref}$ be natural mortality at a reference length $l_{ref}$. The following options are supported:

---

**Constant mortality**

Natural mortality is assumed to be constant across all lengths:

$$
M(l) = M_{ref}
$$

---

**Inverse length scaling**

Natural mortality declines inversely with length:

$$
M(l) = M_{ref} \left(\frac{l_{ref}}{l}\right)
$$

This simple formulation captures the general expectation that smaller individuals experience higher mortality.

---

**Lorenzen (2000)**

Natural mortality scales with body weight $W(l)$, typically approximated from a length–weight relationship:

$$
M(l) = M_{ref} \left(\frac{W(l)}{W(l_{ref})}\right)^{-0.288}
$$

where:

$$
W(l) = a l^b
$$

This formulation implies a smooth decline in mortality with increasing size and is widely used in length-based and ecosystem models.

---

**Gislason et al. (2010)**

Natural mortality is modelled as a function of length and asymptotic length $L_\infty$:

$$
\log M(l) = -0.55 - 1.61 \log l + 1.44 \log L_\infty + \log k
$$



where $l$ is length, $L_\infty$ is asymptotic length, and $k$ is the von Bertalanffy growth coefficient.

This formulation captures empirical scaling relationships across species and implies that natural mortality decreases with size but increases with growth rate and maximum size.


---

These alternative formulations allow sensitivity analyses to assumptions about size-dependent mortality, which can influence estimates of selectivity, SPR, and equilibrium reference points.


## Selectivity

`FLicc` currently supports simple parametric selectivity functions that can also be combined into mixtures:

**Logistic**

$$
sel_i = \frac{1}{1 + \exp(-S_s(L_i - S_m))}
$$

**Normal**

$$
sel_i = \exp\{-S_s(L_i - S_m)^2\}
$$

**Double-sided normal**

$$
sel_i =
\begin{cases}
\exp\{-S_{s1}(L_i - S_m)^2\}, & L_i < S_m \\
\exp\{-S_{s2}(L_i - S_m)^2\}, & L_i > S_m
\end{cases}
$$

These functions correspond roughly to asymptotic, symmetric dome-shaped, and flexible dome-shaped patterns. The ability to fit dome-shaped selectivity is especially important in gillnet and mixed-gear fisheries, where assuming logistic selectivity alone can bias SPR estimates. Medley’s simulations showed much poorer performance when dome-shaped selectivity was mis-specified as logistic. 

A key practical advantage of fitting multiple gears simultaneously is that one gear may still detect larger fish even when another gear is dome-shaped or truncated. That improves the ability to separate selectivity from the declining abundance pattern with length. :contentReference[oaicite:8]{index=8}

---

## Basic population dynamics and equilibrium outputs

The model is equilibrium-based, but `FLicc` extends this into a practical FLR workflow. The package now supports:

- fitting multiple years in a single framework,
- estimating common biological and selectivity parameters across years,
- allowing fishing mortality to vary by year,
- equilibrium reference curve functions for SPR, SSB/SSB0 and relative yield,
- FLR-style extraction of fitted quantities and diagnostics. 

In this sense, `FLicc` keeps the equilibrium interpretation of the original model, but expands it into a multi-year estimation platform that can stabilize inference on shared parameters while still allowing annual variation in exploitation. This is particularly useful when annual data are individually weak, but collectively informative.

---

## What is new relative to fishblicc?

Relative to the original `fishblicc` implementation in Stan, `FLicc` introduces several key practical and methodological extensions:

- **Multi-year estimation framework with TMB**  
  The model is implemented in Template Model Builder (TMB), enabling efficient joint estimation across multiple years. This allows the model to scale to larger datasets and more complex structures while retaining fast optimization and access to automatic differentiation for uncertainty estimation. Joint fitting across years improves parameter identifiability by leveraging temporal replication in the data.

- **Structured variation in fishing mortality**  
  Annual fishing mortality is modelled as a log-scale random walk:
  
  
$$
\log F_y = \log F_{y-1} + \epsilon_y, \quad \epsilon_y \sim \mathcal{N}(0, \sigma_F^2)
$$
  
  
  This introduces temporal structure that stabilises estimation when fitting multiple years jointly, allowing fishing mortality to evolve smoothly while still capturing interannual variability. It also facilitates the estimation of life-history and selectivity parameters by borrowing strength across years.

- **Improved inference through multi-gear data integration**  
  When multiple gears are included, the model exploits differences in selectivity patterns to better resolve fishing mortality and selectivity. In particular, gears that sample different parts of the size or age range—especially those retaining larger individuals—help reduce confounding between mortality and selectivity, improving robustness of parameter estimates.

- **Equilibrium model functions**  
  The package provides equilibrium-based diagnostics and reference point tools, including SPR curves, relative biomass curves, and yield curves. These enable direct integration into management workflows, support scenario testing, and facilitate compatibility with FLR-based pipelines.
---

## Length-based indicators: from large fish to stock structure

A second area of extension in `FLicc` is the development of length-based indicators of stock structure.

### LBIspr concept

The idea behind `LBIspr` is to compare the observed proportion of “large” fish to the proportion expected under a target SPR level, in close analogy to age-based indicators defined relative to $F_{MSY}$. Griffiths et al. defined the age-based indicator $ABI_{MSY}$ as:

$$
ABI_{MSY,t} = \frac{P_t}{P_{MSY}}
$$

where $P_t$ is the proportion of fish above a reference age in year $t$, and $P_{MSY}$ is the corresponding equilibrium proportion above that threshold under $F_{MSY}$. Their results showed that $ABI_{MSY}$ tracked exploitation pressure, was negatively related to $F/F_{MSY}$, and had useful classification skill relative to $B/B_{MSY}$. 

`LBIspr` uses the same logic in a length-based setting. Rather than defining the reference threshold from age structure at $F_{MSY}$, the threshold is defined from the equilibrium length structure under a chosen SPR target, for example SPR40. Let:

- $(L_{SPRx}$ be the length threshold associated with the target SPR level,
- $(P_t(L \ge L_{SPRx})$ be the observed proportion of fish above that threshold,
- $(P_{SPRx}$ be the expected equilibrium proportion above that threshold under the target SPR level.

Then a simple length-based index can be written as:

$$
LBI_{SPRx,t} = \frac{P_t(L \ge L_{SPRx})}{P_{SPRx}}
$$

Values below 1 indicate fewer large fish than expected under the target structure; values above 1 indicate more. This gives an interpretable structure-based indicator that can be tracked through time and related to both exploitation and target SPR conditions.

### Conceptual link to Griffiths et al. and Goodyear

This idea builds directly on the equilibrium-reference concept of $ABI_{MSY}$: observed structure is compared to expected structure under a management reference point. 

It also connects to earlier conceptual work by Goodyear, who proposed NZ50 as a way to monitor the frequency of large individuals in catches. NZ50 is the smallest sample size required for a fish at or above a threshold size to appear in 50% of random samples. He showed that:

$$
NZ50 = \frac{\log(0.5)}{\log(p)}
$$

where $p$ is the cumulative probability at the threshold. Goodyear argued that such threshold-based metrics can track the relative abundance of large fish, are independent of the overall shape of the catch-length distribution, and are sensitive to changes in fishing mortality when catch removes multiple ages in long-lived species. He also suggested that explicit treatment of the largest fish in the catch would be a useful adjunct to standard reference points such as $B/B_{MSY}$ and $F/F_{MSY}$. 

In that sense, `LBIspr` sits naturally between these two lines of work:

- like Goodyear, it focuses on structure in the upper tail of the size distribution, and
- like Griffiths et al., it benchmarks that structure against an equilibrium management target.

---


## Installation

`FLicc` depends on `FLCore`, `ggplotFL`, and `TMB`. The FLR packages are available from the FLR r-universe repository, while `TMB` is installed from CRAN. On Windows, installing `TMB` from source requires Rtools.



### 1. Install required packages

```r
install.packages("devtools")
install.packages("TMB", type = "source")
install.packages("ggplot2")
```

`TMB` is the estimation engine used by `FLicc`. Installing from source ensures compatibility when compiling models and linking against the local toolchain.

### 2. Windows only: install Rtools (if not already installed)

On Windows, `TMB` requires Rtools. You can install in R from:


```r
browseURL("https://cran.r-project.org/bin/windows/Rtools/")
```

After installation, restart R before proceeding.


### 3. Install FLR dependencies

```r
install.packages(
  c("FLCore", "ggplotFL"),
  repos = c("https://flr.r-universe.dev", "https://cloud.r-project.org")
)
```

### 4. Install FLicc

```r
devtools::install_github("henning-winker/FLicc")
```


### 5. Load the package

```r
library(FLicc)
```

## Typical workflow

A standard `FLicc` analysis consists of:

1. Preparing gear-specific length-frequency data and life-history inputs.
2. Building an `FLStockLen` object with the chosen natural mortality model.
3. Checking life-history assumptions such as weight-at-length, maturity-at-length, and mortality-at-length.
4. Fitting one or more years jointly.
5. Inspecting fitted versus observed length compositions and selectivity patterns.
6. Converting fitted outputs to FLR-compatible objects.
7. Evaluating equilibrium reference curves and SPR-based indicators.

For a full reproducible example, see the test script:
[FLicc test workflow](https://github.com/Henning-Winker/FLicc/blob/main/test.R)

---



## Short workflow example

```r


library(FLicc)


data("alfonsino")

# Example input date structure
LFD.df <- lfd_long_to_wide(as.data.frame(lfd_alfonsino))

head(LFD.df)
unique(LFD.df$gear)
unique(LFD.df$len)

# Convert to FLQuants
lfd <- FLQuantLen(LFD.df,unit="cm",midL=FALSE)
# Observed Frequencies
plot_lfd(lfd,type="relmax")

# Normalized to maximum
plot_lfd(lfd)


# Specify Life History
lhpar <- FLPar(
  linf = 55.7,
  k    = 0.08,
  M    = 0.162,
  L50  = 31.1859,
  a    = 0.004721956/1000,
  b    = 3.146168
)


# Build FLStockLen input
stklen <- stocklen(lfd,lhpar,m_model="constant")

# check life history input
plot_lw(stklen)
plot_mat(stklen)
plot_m(stklen)

# Check m model shapes
m_models = c("constant", "inverse", "Lorenzen", "Gislason")
# loop through

m_stks <- FLStocks(lapply(m_models,function(x){
  stocklen(lfd,lhpar,m_model=x)
}))
names(m_stks) <- m_models
plot_m(m_stks)


# Fit model
fit <- fiticc(lfd, stklen,sel_fun=c("dsnormal","logistic"),catch_by_gear =c(0.7,0.3),
              settings=list(prior_sigmaF = c(log(0.5), 0.3,1)))
# log-likelihood
ll <- LLflicc(fit)
ll[[1]]
AIC(ll)
BIC(ll)
# convergence check
flicc_convergence(fit)
# internal optimization to reduce gradient
fit$restart_log

plot_spr(fit)
# FLReport structure
summary(fit$report)
fcur_flicc(fit)
sprcur_flicc(fit)

# Plotting observed vs predicted
plot_len(fit)
plot_len(fit,by_gear = T,year=2020:2024)

# Plot Selectivity by gear
plot_sel(fit)
selpars_flicc(fit)
do.call(c,lapply(selpars_flicc(fit),function(x)x))


# Option to create FLStockLen
stkl <- flicc_stklen(fit)
# Plot, e.g., fishery selectivity weighted by the ratio of catches
plot_sel(stkl)


# Equilibrium dynamics
eqstk <- eqstklen(fit,s=0.75)
eqstk@refpts
plot_eqcurves(eqstk)

# get repts
## compute per-recruit brp
prbrp <- prbrp_flicc(fit)
## compute Fmsy as function of steepness
fmsy_flicc(prbrp,s=0.75)
## compare to Fspr40 and Fspr35
fspr_flicc(fit,spr=40)
fspr_flicc(fit,spr=35)
## control
round(spr_flicc(fit,F=fspr_flicc(fit,spr=40)),3)


# LBIspr indicator
LBIidx <- LBIspr(fit,thresh = 0.75)

# Plot LBIspr with FLicc
plot_LBIspr(fit,thresh = 0.75)


# Status
plot_spr(fit)
plot_lbfao(fit)

# Convert to simplified FLStockR
stk <- flicc2FLStockR(fit)
stk@refpts
plot_LBAdvice(stk)
# relative to MSY proxy
stkr <- flicc2FLStockR(fit,rel=T)
plot_LBAdvice(stkr,panel=1)+ylim(0.,1.5)


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Fit model each year separately (LBSPR-like)
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

fit.y <- fiticc(lfd, stklen,sel_fun=c("dsnormal","logistic"),catch_by_gear =c(0.7,0.3),
              settings=list(prior_sigmaF = c(log(0.5), 0.3,1)),by_year=TRUE)

stky <- flicc_stklen(fit.y)
# Plot fishery selectivity estimated for each year
plot_sel(stky)

# compare
plot_spr(list(all.yr=fit,each.y=fit.y))

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Compare fishblicc and LBSPR population and error models
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
#> LBSRP: pop_model = "gtg", obs_model = "mn" # multinomial (Default)
#> fishblicc: pop_model = "gamma", obs_model = "nb" # negative bionomial

#> Note ESS is based on the number of observations (like in LBSRP,fishblicc).
#> This can be adjusted, for example, by:
#>
lfd.ess <- lfdess(lfd,ess.g=c(Trawl=250,Gillnet=350))
plot_lfd(lfd.ess)

# Slower
system.time({
fit.gamma.nb <- fiticc(lfd.ess, stklen,sel_fun=c("dsnormal","logistic"),catch_by_gear =c(0.7,0.3),
                settings=list(pop_model="gamma",obs_model="nb"))
})

# Faster
system.time({
  fit.gtg.mn <- fiticc(lfd.ess, stklen,sel_fun=c("dsnormal","logistic"),catch_by_gear =c(0.7,0.3),
                   settings=list(pop_model="gtg",obs_model="mn"))
})

# compare
plot_spr(list(gamma.nb=fit.gamma.nb,
              gtg.mn=fit.gtg.mn))


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# fishblicc test example
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

data("fishblicc_example")

# Build stock with inverse M model
stklen <- stocklen(lfd_fishblicc,lhpar_fishblicc,m_model="inverse")

plot_m(stklen)

# Fit model with penalties (priors) on Linf, MK and CVL approximating fishblicc
# pop_model = "gamma"
# obs_model="nb"
# linf.sd=2/40   (fishblicc sdev = 2, here as log.se ~ sdev/linf)
# Mk.sd=0.1   (also lognormal in fishblicc)
# CVL.sd=0.1   (parameterized as Galpha in fishblicc)


fit <- fiticc(lfd_fishblicc, stklen,sel_fun=c("dsnormal","dsnormal","dsnormal"),catch_by_gear =c(0.1802070, 0.2101353, 0.6096577),
              settings=list(pop_model = "gamma", obs_model="nb",CVL=0.14,CVL.sd=0.1,linf.sd=2/40,Mk.sd=0.1))

# Convergence
flicc_convergence(fit)
sprcur_flicc(fit)
# flicc spr - 0.3408
# fishblicc spr = 0.340

plot_sel(fit)

# Plotting observed vs predicted
plot_len(fit,by_gear = T)
plot_len(fit)


# Equilibrium dynamics
eqstk <- eqstklen(fit,s=0.7)
eqstk@refpts
plot_eqcurves(eqstk)



```




## Future developments

Future developments could include:

- exploring random effects on selectivity parameters, for example to allow selected gear-specific selectivity parameters to vary among years while estimating shared mean patterns,
- isopleth plots for visualising equilibrium trade-offs among fishing mortality, SPR, biomass and yield,
- additional selectivity options 
- diagnostics: residuals, likelihood profiling, sensitivities
- MSE implementation: develop management procedure mp.flicc with harvest control rule options 

---

## References

- Griffiths, C.A., Winker, H., Bartolino, V., Wennhage, H., Orio, A. and Cardinale, M., 2024. Including older fish in fisheries management: A new age‐based indicator and reference point for exploited fish stocks. Fish and Fisheries, 25(1), pp.18-37.
- Goodyear, C. P. 2015a. Understanding maximum size in the catch: Atlantic blue marlin as an example. Transactions of the American Fisheries Society, 144:274-282.
- Hordyk, A., Ono, K., Valencia, S., Loneragan, N. and Prince, J., 2015. A novel length-based empirical estimation method of spawning potential ratio (SPR), and tests of its performance, for small-scale, data-poor fisheries. ICES Journal of Marine Science, 72(1), pp.217-231.
- Hordyk, A.R., Ono, K., Prince, J.D. and Walters, C.J., 2016. A simple length-structured model based on life history ratios and incorporating size-dependent selectivity: application to spawning potential ratios for data-poor stocks. Canadian Journal of Fisheries and Aquatic Sciences, 73(12), pp.1787-1799.
- Kell, L.T., Mosqueira, I., Grosjean, P., Fromentin, J.M., Garcia, D., Hillary, R., Jardim, E., Mardle, S., Pastoors, M.A., Poos, J.J. and Scott, F., 2007. FLR: an open-source framework for the evaluation and development of management strategies. ICES Journal of Marine Science, 64(4), pp.640-646.
- Kristensen, K., Nielsen, A., Berg, C.W., Skaug, H. and Bell, B.M., 2016. TMB: automatic differentiation and Laplace approximation. Journal of statistical software, 70, pp.1-21.
- Medley, P.A., 2025. A new Bayesian catch curve stock assessment model for the analysis of length data from multi-gear fisheries. ICES Journal of Marine Science, 82(12), p.fsaf224.
