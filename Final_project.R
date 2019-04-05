library("ReacTran")
library("marelac")
if (!require("viridis")) {
    install.packages("viridis")
    library("viridis")
}
if (!require("lubridate")) {
    install.packages("lubridate")
    library("lubridate")
}

start_time <- Sys.time()
# setwd("Final project/")

nitrate_rhine <- read.csv2('data/NO3_CONCTTE_LOBPTN.csv')
nitrate_rhine$datum_tijd <- dmy_hms(paste(nitrate_rhine$Datum,
                                          nitrate_rhine$Tijd))
nitrate_rhine$hours <- as.double(nitrate_rhine$datum_tijd -
                                      dmy_hm('1-1-2010 0:00'),
                                  units='hours')
nitrate_rhine$meting_umol <- nitrate_rhine$Meting / 14 * 1000

ammon_rhine <- read.csv2('data/NH4_CONCTTE_LOBPTN.csv')
ammon_rhine$datum_tijd <- dmy_hms(paste(ammon_rhine$Datum,
                                          ammon_rhine$Tijd))
ammon_rhine$hours <- as.double(ammon_rhine$datum_tijd -
                                     dmy_hm('1-1-2010 0:00'),
                                 units='hours')
ammon_rhine$meting_umol <- ammon_rhine$Meting / 14 * 1000

DOM_rhine <- read.csv2('data/Corg_CONCTTE_LOBPTN.csv')
DOM_rhine$datum_tijd <- dmy_hms(paste(DOM_rhine$Datum,
                                        DOM_rhine$Tijd))
DOM_rhine$hours <- as.double(DOM_rhine$datum_tijd -
                                   dmy_hm('1-1-2010 0:00'),
                               units='hours')
DOM_rhine$meting_umol <- DOM_rhine$Meting / 12 / 106 * 1000

temp_rhine <- read.csv2('data/20190403_034.csv')
temp_rhine$datum_tijd <- dmy_hms(paste(temp_rhine$WAARNEMINGDATUM,
                                       temp_rhine$WAARNEMINGTIJD))
temp_rhine$hours <- as.double(temp_rhine$datum_tijd-
                                  dmy_hm('1-1-2010 0:00'),
                              units='hours')
temp_rhine$hours[1] <- 0
temp_rhine$temp <- temp_rhine$NUMERIEKEWAARDE
temp_rhine <- temp_rhine[,c('hours', 'temp')]

givendata <- read.table(
    file = "data/p4-data.dat",
    sep = "\t",
    dec = ".",
    col.names = c("distance", "O2", "NH3")
)
distance <- givendata[[1]]
dissolvedO2 <- givendata[[2]]
amonium <- givendata[[3]]


## Affinity constants
kO2 <- 20 # [micromol/L]
kNO3 <- 35 # [micromol/L]

# parameters
L <- 500 # [m]
N <- 500
por <- 0.4 # [m^3_pw/m^3]
# tort <- 1 - log(por^2)
u <- 0.1 # [m/h]
dispersivity <- 1.5 # [m]
D <- u * dispersivity # [m^2/h]

# grids
grid <- setup.grid.1D(x.up = 0, N = N, L = L)
por.grid <- setup.prop.1D(value = por, grid = grid)
D.grid <- setup.prop.1D(value = D, grid = grid)
u.grid <- setup.prop.1D(value = u, grid = grid)
svf <- 1 - por
svf.grid <- setup.prop.1D(value = svf, grid = grid)


The_model <- function(t, state, parms) {
    with(as.list(c(parms)), {
        sat_conc <- gas_satconc(S = S, t = Temp, P = 1, species = c("O2", "N2"))
        SN2 <- as.numeric(sat_conc["N2"]) # [micromol/L]
        SO2 <- as.numeric(sat_conc["O2"]) # [micromol/L]

                # Initialisation of state variables
        DOM <- state[1:N]
        O2 <- state[(N + 1):(2 * N)]
        NH3 <- state[(N * 2 + 1):(3 * N)]
        NO3 <- state[(N * 3 + 1):(4 * N)]
        N2 <- state[(N * 4 + 1):(5 * N)]

        # Transport
        tranDOM <- tran.1D(
            C = DOM, C.up = DOMin,
            D = D.grid, v = u.grid,
            VF = por.grid, dx = grid
        )
        tranO2 <- tran.1D(
            C = O2, C.up = O2in,
            D = D.grid, v = u.grid,
            VF = por.grid, dx = grid
        )
        tranNH3 <- tran.1D(
            C = NH3, C.up = NH3in,
            D = D.grid, v = u.grid,
            VF = por.grid, dx = grid
        )
        tranNO3 <- tran.1D(
            C = NO3, C.up = NO3in,
            D = D.grid, v = u.grid,
            VF = por.grid, dx = grid
        )
        tranN2 <- tran.1D(
            C = N2, C.up = SN2,
            D = D.grid, v = u.grid,
            VF = por.grid, dx = grid
        )

        # Reactions
        Rorg.rem <- k1 * DOM
        Raero <- Rorg.rem * (O2 / (O2 + kO2))
        Rdnit <- Rorg.rem * (NO3 / (kNO3 + NO3)) * (kO2 / (kO2 + O2))
        Rnit <- k2 * NH3 * O2
        Raeration <- (k3 * SO2) * (1 - O2 / SO2)
        Rdgas <- k3 * SN2 * (1 - (N2 / SN2))

        O2_per_DOM <- 106
        NH3_per_DOM <- 16
        NO3_per_DOM <- 424 / 5
        NO3_per_NH3 <- 1
        N2_per_DOM <- 212 / 5
        O2_per_NH3 <- 2

        # total reactions per component
        reacDOM <- -Raero - Rdnit
        reacO2 <- -Raero * O2_per_DOM - Rnit * O2_per_NH3 + Raeration
        reacNH3 <- +Raero * NH3_per_DOM + Rdnit * NH3_per_DOM - Rnit
        reacNO3 <- -Rdnit * NO3_per_DOM + Rnit * NO3_per_NH3
        reacN2 <- +Rdnit * N2_per_DOM + Rdgas

        dCdtDOM <- tranDOM$dC + reacDOM
        dCdtO2 <- tranO2$dC + reacO2
        dCdtNH3 <- tranNH3$dC + reacNH3
        dCdtNO3 <- tranNO3$dC + reacNO3
        dCdtN2 <- tranN2$dC + reacN2

        # Assemble the total rate of change
        return(list(
            c(
                dCdtDOM = dCdtDOM,
                dCdtO2 = dCdtO2,
                dCdtNH3 = dCdtNH3,
                dCdtNO3 = dCdtNO3,
                dCdtN2 = dCdtN2
            ),
            Raero = Raero,
            Rdenit = Rdnit,
            Rnit = Rnit,
            Raeration = Raeration,
            Rdgas = Rdgas
        ))
    })
}

trans_model <- function(t, state, parms) {
    with(as.list(c(parms)), {

        Temp <- approx(temp_rhine$hours, temp_rhine$temp, t)$y
        sat_conc <- gas_satconc(S = S, t = Temp, P = 1, species = c("O2", "N2"))
        SN2 <- as.numeric(sat_conc["N2"]) # [micromol/L]
        SO2 <- as.numeric(sat_conc["O2"]) # [micromol/L]

                # Initialisation of state variables
        DOM <- state[1:N]
        O2 <- state[(N + 1):(2 * N)]
        NH3 <- state[(N * 2 + 1):(3 * N)]
        NO3 <- state[(N * 3 + 1):(4 * N)]
        N2 <- state[(N * 4 + 1):(5 * N)]

        # Transport
        tranDOM <- tran.1D(
            C = DOM,
            C.up = approx(DOM_rhine$hours, DOM_rhine$meting_umol, t)$y, # - cos((t)/365/24*2*pi) + 2e-4*t*(t < (365*24*3)),
            D = D.grid, v = u.grid,
            VF = por.grid, dx = grid
        )
        tranO2 <- tran.1D(
            C = O2, C.up = O2in,
            D = D.grid, v = u.grid,
            VF = por.grid, dx = grid
        )
        tranNH3 <- tran.1D(
            C = NH3,
            C.up = approx(ammon_rhine$hours, ammon_rhine$meting_umol, t)$y,
            D = D.grid, v = u.grid,
            VF = por.grid, dx = grid
        )
        tranNO3 <- tran.1D(
            C = NO3,
            # C.up = NO3in,
            C.up = approx(nitrate_rhine$hours, nitrate_rhine$meting_umol, t)$y,
            D = D.grid, v = u.grid,
            VF = por.grid, dx = grid
        )
        tranN2 <- tran.1D(
            C = N2, C.up = SN2,
            D = D.grid, v = u.grid,
            VF = por.grid, dx = grid
        )

        # Reactions
        Rorg.rem <- k1 * DOM
        Raero <- Rorg.rem * (O2 / (O2 + kO2))
        Rdnit <- Rorg.rem * (NO3 / (kNO3 + NO3)) * (kO2 / (kO2 + O2))
        Rnit <- k2 * NH3 * O2
        Raeration <- (k3 * SO2) * (1 - O2 / SO2)
        Rdgas <- k3 * SO2 * (1 - (N2 / SN2))

        O2_per_DOM <- 106
        NH3_per_DOM <- 16
        NO3_per_DOM <- 424 / 5
        NO3_per_NH3 <- 1
        N2_per_DOM <- 212 / 5
        O2_per_NH3 <- 2

        # total reactions per component
        reacDOM <- -Raero - Rdnit
        reacO2 <- -Raero * O2_per_DOM - Rnit * O2_per_NH3 + Raeration
        reacNH3 <- +Raero * NH3_per_DOM + Rdnit * NH3_per_DOM - Rnit
        reacNO3 <- -Rdnit * NO3_per_DOM + Rnit * NO3_per_NH3
        reacN2 <- +Rdnit * N2_per_DOM + Rdgas

        dCdtDOM <- tranDOM$dC + reacDOM
        dCdtO2 <- tranO2$dC + reacO2
        dCdtNH3 <- tranNH3$dC + reacNH3
        dCdtNO3 <- tranNO3$dC + reacNO3
        dCdtN2 <- tranN2$dC + reacN2

        # Assemble the total rate of change
        return(list(
            c(
                dCdtDOM = dCdtDOM,
                dCdtO2 = dCdtO2,
                dCdtNH3 = dCdtNH3,
                dCdtNO3 = dCdtNO3,
                dCdtN2 = dCdtN2
            ),
            Raero = Raero,
            Rdenit = Rdnit,
            Rnit = Rnit,
            Raeration = Raeration,
            Rdgas = Rdgas
        ))
    })
}

solve_dom_model <- function(O2in = 210, DOMin = 4.71,
                            NH3in = 0, NO3in = 100,
                            k1 = 1e-3, k2 = 5e-4, k3 = 5e-4,
                            S = 0, t = 10) {

    # parameters
    parms <- c(
        O2in = O2in,
        DOMin = DOMin,
        NH3in = NH3in,
        NO3in = NO3in,
        k1 = k1,
        k2 = k2,
        k3 = k3,
        kO2 = kO2,
        kNO3 = kNO3,
        S = S,
        Temp = t
    )

    state <- rep(0, length.out = (N * 5))

    std <- steady.1D(
        y = state,
        func = The_model,
        parms = parms,
        names = c("DOM", "O2", "NH3", "NO3", "N2"),
        nspec = 5,
        pos = TRUE
    )
    return(std)
}

solve_dom_model_transient <- function(times,
                                      O2in = 210, DOMin = 4.71,
                                      NH3in = 0, NO3in = 100,
                                      k1 = 1e-3, k2 = 5e-4, k3 = 5e-4,
                                      S = 0, t = 10) {

    # parameters
    parms <- c(
        O2in = O2in,
        DOMin = DOMin,
        NH3in = NH3in,
        NO3in = NO3in,
        k1 = k1,
        k2 = k2,
        k3 = k3,
        kO2 = kO2,
        kNO3 = kNO3,
        S = S,
        Temp = t
    )

    state <- rep(0, length.out = (N * 5))

    std <- steady.1D(
        y = state,
        func = The_model,
        parms = parms,
        names = c("DOM", "O2", "NH3", "NO3", "N2"),
        nspec = 5,
        pos = TRUE
    )
    state <- std$y

    out <- ode.1D(
        y = state,
        times = times,
        func = trans_model,
        parms = parms,
        names = c("DOM", "O2", "NH3", "NO3", "N2"),
        nspec = 5
    )
    return(out)
}

# plotting model
std <- solve_dom_model(k1 = 2e-3, k2 = 3.5e-4, k3 = 3e-4)
obs <- matrix(data = NA, nrow = nrow(givendata), ncol = 6, byrow = FALSE)
colnames(obs) <- c("distance", "DOM", "O2", "NH3", "NO3", "N2")
obs[, "distance"] <- distance
obs[, "O2"] <- dissolvedO2
obs[, "NH3"] <- amonium
plot(std, obs = obs, grid = grid$x.mid)

rates <- cbind(std$Raero, std$Rdenit, std$Rnit, std$Raeration, std$Rdgas)
matplot(grid$x.mid, rates, type = 'l')
legend('topright',
       c('Aero', 'Denit', 'Nit', 'Aera', 'Degas'),
       lty = 1:5,
       col = 1:5)


max_DOM <- 3 / 12 / 106 * 1000
max_NO3 <- 25 / (14 + 3 * 16) * 1000
max_NH3 <- 0.05 / (14 + 3 * 1) * 1000

in_bounds <- (
    std$y[, "DOM"] < max_DOM &
        std$y[, "NO3"] < max_NO3 &
        std$y[, "NH3"] < max_NH3
)
first_in_bounds <- min(grid$x.mid[in_bounds])

total_time = 24*365*8 # 8 years in hours
step_size = 24*14
times <- seq(0, total_time, by = step_size)
out <- solve_dom_model_transient(times, k1 = 2e-3, k2 = 3.5e-4, k3 = 3e-4)
out[,'time'] <- out[,'time']/24/365
image(out, legend = TRUE, grid = grid$x.mid, col = viridis(256), xlab = 'Time [yr]')
image(out, method='persp', grid = grid$x.mid)#, legend = TRUE, col = viridis(256), xlab = 'Time [yr]')
print(Sys.time() - start_time)
# out_diff <- out
# out_diff[,2:2501] <- out_diff[,2:2501] - rep(std$y, each = 1826)
# image(out_diff, legend = TRUE, grid = grid$x.mid, col = RColorBrewer::brewer.pal(11, 'Spectral'), zlim = c(-150, 150))
