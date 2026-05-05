
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
system.time({
fit <- fiticc(lfd, stklen,sel_fun=c("dsnormal","logistic"),catch_by_gear =c(0.7,0.3),
              settings=list(prior_sigmaF = c(log(0.5), 0.3,1),CVL=0.1),n_restart = 10)
})
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

lhpar_fishblicc <- rbind(lhpar_fishblicc,
          FLPar(M=lhpar_fishblicc["k"]*lhpar_fishblicc["Mk"]))
# Build stock with inverse M model
stklen <- stocklen(lfd_fishblicc,lhpar_fishblicc,m_model="inverse",reflen = 21.9)

# Fit model with penalties (priors) on Linf, MK and CVL approximating fishblicc
# pop_model = "gamma"
# obs_model="nb"
# linf.sd=2/40   (fishblicc sdev = 2, here as log.se ~ sdev/linf)
# Mk.sd=0.1   (also lognormal in fishblicc)
# CVL.sd=0.1   (parameterized as Galpha in fishblicc)


fit <- fiticc(lfd_fishblicc, stklen,sel_fun=c("dsnormal","dsnormal","dsnormal"),catch_by_gear =c(0.1802070, 0.2101353, 0.6096577),
              settings=list(pop_model = "gamma", obs_model="nb",
                            CVL=0.14,
                            CVL.sd=0.1,
                            linf.sd=2/40,
                            Mk.sd=0.1))

stkl <- flicc_stklen(fit)
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

