# Equilibrium FLicc on OM-like simulations: lag and sensitivity

Alfonsino-like OM (fixed-CV length-at-age, annual pulse recruitment, BH h = 0.75, sigmaR 0.5), two gears,
two-stage sampling (n_age 50/30, n_len 150/100), length data 2012-2024, 10 iterations. Tail cuts are absolute lengths
(50.1/47.3/44.5/41.7 cm = 0.90/0.85/0.80/0.75 x true Linf). Relative error in SPR at the year's F.

## A. Lag under changing F (correct biology, logistic Gillnet)

slope_ratio = estimated / true linear trend in SPR over the data years (1 = no lag; rebuild and decline only).
noise = SD of the year-to-year change in the error.

| scen | model | cut | conv | bias | rmse | bias_last3 | slope_ratio | noise |
|---|---|---|---|---|---|---|---|---|
| decline | sF_est | none | 10/10 | +217% | +306% | +275% |  1.54 | 0.056 |
| decline | sF_est | c0.90 | 10/10 | +245% | +330% | +282% |  1.49 | 0.057 |
| decline | sF_est | c0.85 | 10/10 | +222% | +311% | +250% |  1.76 | 0.057 |
| decline | sF_est | c0.80 | 10/10 | +194% | +262% | +193% |  1.87 | 0.054 |
| decline | sF_est | c0.75 | 10/10 | +135% | +202% | +125% |  1.56 | 0.054 |
| decline | sF0.3 | none | 10/10 | +238% | +317% | +271% |  2.37 | 0.059 |
| decline | sF0.3 | c0.90 | 10/10 | +260% | +336% | +269% |  2.46 | 0.061 |
| decline | sF0.3 | c0.85 | 10/10 | +240% | +318% | +217% |  2.56 | 0.060 |
| decline | sF0.3 | c0.80 | 10/10 | +205% | +270% | +168% |  2.55 | 0.060 |
| decline | sF0.3 | c0.75 | 10/10 | +130% | +209% | +111% |  2.30 | 0.059 |
| decline | sF0.5 | none | 10/10 | +252% | +336% | +254% |  2.75 | 0.079 |
| decline | sF0.5 | c0.90 | 10/10 | +254% | +356% | +248% |  2.99 | 0.078 |
| decline | sF0.5 | c0.85 | 10/10 | +234% | +339% | +196% |  3.10 | 0.078 |
| decline | sF0.5 | c0.80 | 10/10 | +208% | +289% | +135% |  3.10 | 0.076 |
| decline | sF0.5 | c0.75 | 10/10 | +142% | +232% | +99% |  2.61 | 0.074 |
| decline | sF1 | none | 10/10 | +239% | +357% | +224% |  3.01 | 0.109 |
| decline | sF1 | c0.90 | 10/10 | +262% | +376% | +236% |  3.31 | 0.117 |
| decline | sF1 | c0.85 | 10/10 | +237% | +359% | +190% |  3.45 | 0.116 |
| decline | sF1 | c0.80 | 10/10 | +217% | +311% | +136% |  3.45 | 0.109 |
| decline | sF1 | c0.75 | 10/10 | +151% | +260% | +93% |  2.94 | 0.109 |
| decline | noRW | none | 10/10 | +256% | +403% | +219% |  3.33 | 0.223 |
| decline | noRW | c0.90 | 10/10 | +273% | +421% | +238% |  3.45 | 0.216 |
| decline | noRW | c0.85 | 10/10 | +269% | +408% | +194% |  3.61 | 0.215 |
| decline | noRW | c0.80 | 10/10 | +241% | +376% | +122% |  3.76 | 0.220 |
| decline | noRW | c0.75 | 10/10 | +174% | +350% | +93% |  3.37 | 0.237 |
| flat | sF_est | none | 10/10 | +7% | +31% | +19% | -0.22 | 0.063 |
| flat | sF_est | c0.90 | 10/10 | +20% | +42% | +28% | -0.25 | 0.063 |
| flat | sF_est | c0.85 | 10/10 | +22% | +46% | +33% | -0.27 | 0.063 |
| flat | sF_est | c0.80 | 10/10 | +27% | +46% | +31% | -0.09 | 0.063 |
| flat | sF_est | c0.75 | 10/10 | +15% | +41% | +22% |  0.27 | 0.063 |
| flat | sF0.3 | none | 10/10 | +5% | +39% | +15% | -0.33 | 0.064 |
| flat | sF0.3 | c0.90 | 10/10 | +19% | +50% | +24% | -0.48 | 0.065 |
| flat | sF0.3 | c0.85 | 10/10 | +23% | +55% | +26% | -0.53 | 0.066 |
| flat | sF0.3 | c0.80 | 10/10 | +24% | +56% | +27% | -0.40 | 0.067 |
| flat | sF0.3 | c0.75 | 10/10 | +12% | +46% | +16% |  0.34 | 0.066 |
| flat | sF0.5 | none | 10/10 | +7% | +43% | +10% | -0.32 | 0.071 |
| flat | sF0.5 | c0.90 | 10/10 | +17% | +56% | +22% | -0.46 | 0.074 |
| flat | sF0.5 | c0.85 | 10/10 | +22% | +62% | +24% | -0.51 | 0.075 |
| flat | sF0.5 | c0.80 | 10/10 | +25% | +64% | +28% | -0.44 | 0.075 |
| flat | sF0.5 | c0.75 | 10/10 | +17% | +54% | +18% |  0.37 | 0.073 |
| flat | sF1 | none | 10/10 | +3% | +48% | +9% | -0.31 | 0.089 |
| flat | sF1 | c0.90 | 10/10 | +16% | +66% | +22% | -0.44 | 0.098 |
| flat | sF1 | c0.85 | 10/10 | +26% | +74% | +24% | -0.49 | 0.101 |
| flat | sF1 | c0.80 | 10/10 | +31% | +78% | +34% | -0.44 | 0.101 |
| flat | sF1 | c0.75 | 10/10 | +20% | +71% | +22% |  0.47 | 0.094 |
| flat | noRW | none | 10/10 | +4% | +59% | +11% | -0.32 | 0.109 |
| flat | noRW | c0.90 | 10/10 | +23% | +95% | +26% | -0.44 | 0.131 |
| flat | noRW | c0.85 | 10/10 | +30% | +107% | +31% | -0.50 | 0.144 |
| flat | noRW | c0.80 | 10/10 | +36% | +118% | +35% | -0.63 | 0.159 |
| flat | noRW | c0.75 | 10/10 | +31% | +115% | +34% |  0.57 | 0.162 |
| rebuild | sF_est | none | 10/10 | -56% | +55% | -54% |  0.14 | 0.079 |
| rebuild | sF_est | c0.90 | 10/10 | -55% | +53% | -53% |  0.18 | 0.079 |
| rebuild | sF_est | c0.85 | 10/10 | -55% | +53% | -54% |  0.21 | 0.079 |
| rebuild | sF_est | c0.80 | 10/10 | -52% | +51% | -51% |  0.26 | 0.079 |
| rebuild | sF_est | c0.75 | 10/10 | -53% | +50% | -50% |  0.19 | 0.078 |
| rebuild | sF0.3 | none | 10/10 | -54% | +55% | -51% |  0.27 | 0.079 |
| rebuild | sF0.3 | c0.90 | 10/10 | -52% | +53% | -48% |  0.35 | 0.080 |
| rebuild | sF0.3 | c0.85 | 10/10 | -53% | +53% | -49% |  0.40 | 0.078 |
| rebuild | sF0.3 | c0.80 | 10/10 | -50% | +51% | -43% |  0.50 | 0.080 |
| rebuild | sF0.3 | c0.75 | 10/10 | -50% | +51% | -46% |  0.47 | 0.079 |
| rebuild | sF0.5 | none | 10/10 | -54% | +55% | -51% |  0.28 | 0.082 |
| rebuild | sF0.5 | c0.90 | 10/10 | -52% | +53% | -47% |  0.37 | 0.085 |
| rebuild | sF0.5 | c0.85 | 10/10 | -52% | +53% | -47% |  0.42 | 0.083 |
| rebuild | sF0.5 | c0.80 | 10/10 | -48% | +51% | -44% |  0.54 | 0.084 |
| rebuild | sF0.5 | c0.75 | 10/10 | -49% | +50% | -42% |  0.51 | 0.085 |
| rebuild | sF1 | none | 10/10 | -54% | +54% | -51% |  0.29 | 0.090 |
| rebuild | sF1 | c0.90 | 10/10 | -51% | +52% | -47% |  0.38 | 0.095 |
| rebuild | sF1 | c0.85 | 10/10 | -50% | +52% | -48% |  0.44 | 0.095 |
| rebuild | sF1 | c0.80 | 10/10 | -48% | +50% | -43% |  0.57 | 0.095 |
| rebuild | sF1 | c0.75 | 10/10 | -47% | +48% | -41% |  0.55 | 0.102 |
| rebuild | noRW | none | 10/10 | -54% | +54% | -52% |  0.30 | 0.096 |
| rebuild | noRW | c0.90 | 10/10 | -50% | +52% | -48% |  0.39 | 0.104 |
| rebuild | noRW | c0.85 | 10/10 | -51% | +52% | -48% |  0.44 | 0.107 |
| rebuild | noRW | c0.80 | 10/10 | -47% | +49% | -40% |  0.59 | 0.120 |
| rebuild | noRW | c0.75 | 10/10 | -45% | +50% | -38% |  0.59 | 0.123 |

## B. Sensitivity: Gillnet dome (fitted as logistic), input Linf and CVL, x tail cut (sigmaF 0.3)

Median bias by assumed Linf multiplier and CVL (rows) and tail cut (columns).

| scen | assumed | bias.none | bias.c0.90 | bias.c0.85 | bias.c0.80 |
|---|---|---|---|---|---|
| flat_dome10 | Linf0.9_CVL0.07 | +72% | +74% | +77% | +72% |
| flat_dome10 | Linf0.9_CVL0.10 | +48% | +57% | +64% | +62% |
| flat_dome10 | Linf0.9_CVL0.13 | +24% | +41% | +48% | +50% |
| flat_dome10 | Linf1.0_CVL0.07 | -11% | -7% | +2% | +5% |
| flat_dome10 | Linf1.0_CVL0.10 | -17% | -12% | -3% | +1% |
| flat_dome10 | Linf1.0_CVL0.13 | -24% | -18% | -8% | -3% |
| flat_dome10 | Linf1.1_CVL0.07 | -44% | -42% | -36% | -31% |
| flat_dome10 | Linf1.1_CVL0.10 | -47% | -44% | -38% | -33% |
| flat_dome10 | Linf1.1_CVL0.13 | -50% | -47% | -40% | -35% |
| flat_dome20 | Linf0.9_CVL0.07 | +114% | +107% | +93% | +88% |
| flat_dome20 | Linf0.9_CVL0.10 | +83% | +85% | +78% | +77% |
| flat_dome20 | Linf0.9_CVL0.13 | +54% | +63% | +61% | +64% |
| flat_dome20 | Linf1.0_CVL0.07 | +4% | +9% | +7% | +14% |
| flat_dome20 | Linf1.0_CVL0.10 | -6% | +3% | +3% | +12% |
| flat_dome20 | Linf1.0_CVL0.13 | -15% | -4% | -3% | +9% |
| flat_dome20 | Linf1.1_CVL0.07 | -39% | -33% | -32% | -25% |
| flat_dome20 | Linf1.1_CVL0.10 | -42% | -35% | -34% | -26% |
| flat_dome20 | Linf1.1_CVL0.13 | -46% | -38% | -37% | -27% |
| flat_logistic | Linf0.9_CVL0.07 | +121% | +110% | +100% | +84% |
| flat_logistic | Linf0.9_CVL0.10 | +91% | +90% | +87% | +74% |
| flat_logistic | Linf0.9_CVL0.13 | +62% | +70% | +72% | +63% |
| flat_logistic | Linf1.0_CVL0.07 | +9% | +12% | +12% | +13% |
| flat_logistic | Linf1.0_CVL0.10 | -1% | +7% | +7% | +9% |
| flat_logistic | Linf1.0_CVL0.13 | -11% | +1% | +2% | +4% |
| flat_logistic | Linf1.1_CVL0.07 | -36% | -30% | -29% | -27% |
| flat_logistic | Linf1.1_CVL0.10 | -40% | -32% | -30% | -28% |
| flat_logistic | Linf1.1_CVL0.13 | -43% | -35% | -32% | -31% |
| rebuild_dome10 | Linf0.9_CVL0.07 | -29% | -29% | -28% | -30% |
| rebuild_dome10 | Linf0.9_CVL0.10 | -36% | -36% | -34% | -34% |
| rebuild_dome10 | Linf0.9_CVL0.13 | -42% | -42% | -40% | -39% |
| rebuild_dome10 | Linf1.0_CVL0.07 | -59% | -58% | -57% | -57% |
| rebuild_dome10 | Linf1.0_CVL0.10 | -61% | -60% | -59% | -59% |
| rebuild_dome10 | Linf1.0_CVL0.13 | -64% | -63% | -62% | -61% |
| rebuild_dome10 | Linf1.1_CVL0.07 | -74% | -73% | -72% | -72% |
| rebuild_dome10 | Linf1.1_CVL0.10 | -75% | -74% | -74% | -73% |
| rebuild_dome10 | Linf1.1_CVL0.13 | -76% | -75% | -75% | -74% |
| rebuild_dome20 | Linf0.9_CVL0.07 | -19% | -21% | -22% | -22% |
| rebuild_dome20 | Linf0.9_CVL0.10 | -30% | -29% | -28% | -27% |
| rebuild_dome20 | Linf0.9_CVL0.13 | -40% | -37% | -35% | -32% |
| rebuild_dome20 | Linf1.0_CVL0.07 | -56% | -55% | -54% | -53% |
| rebuild_dome20 | Linf1.0_CVL0.10 | -59% | -58% | -56% | -54% |
| rebuild_dome20 | Linf1.0_CVL0.13 | -62% | -61% | -59% | -57% |
| rebuild_dome20 | Linf1.1_CVL0.07 | -73% | -72% | -71% | -70% |
| rebuild_dome20 | Linf1.1_CVL0.10 | -71% | -73% | -71% | -71% |
| rebuild_dome20 | Linf1.1_CVL0.13 | -73% | -74% | -73% | -72% |
| rebuild_logistic | Linf0.9_CVL0.07 | -10% | -13% | -14% | -18% |
| rebuild_logistic | Linf0.9_CVL0.10 | -21% | -20% | -20% | -22% |
| rebuild_logistic | Linf0.9_CVL0.13 | -30% | -28% | -27% | -27% |
| rebuild_logistic | Linf1.0_CVL0.07 | -52% | -51% | -49% | -50% |
| rebuild_logistic | Linf1.0_CVL0.10 | -55% | -54% | -51% | -51% |
| rebuild_logistic | Linf1.0_CVL0.13 | -59% | -56% | -54% | -53% |
| rebuild_logistic | Linf1.1_CVL0.07 | -71% | -69% | -68% | -67% |
| rebuild_logistic | Linf1.1_CVL0.10 | -71% | -70% | -69% | -68% |
| rebuild_logistic | Linf1.1_CVL0.13 | -71% | -72% | -70% | -69% |
