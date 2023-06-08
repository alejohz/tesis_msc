#Original
#values <- c(27, 27, 7, 24, 39, 40, 24, 45, 36, 37, 31, 47, 16, 24, 6, 21,
#35, 36, 21, 40, 32, 33, 27, 42, 14, 21, 5, 19, 31, 32, 19, 36,
#29, 29, 24, 42, 15, 24, 21)
# y <- c(27, 27, 7, 24, 39, 40, 24, 45, 36, 37, 31, 47, 16, 24,
# 6, 21, 35, 36, 21, 40, 32, 33, 27, 42, 14, 21, 5, 19, 31, 32, 19,
# 36, 29, 29, 24, 42, 15, 24, 21, 35)
library(robets)
library(dplyr)

outlier_multi <- function(x) {
    if (x > 0) {
        return(1 + x)
    } else {
        return(exp(x))
    }
}

model <- function(type, out = FALSE) {
    ts <- 40
    m <- 4
    sigma <- 0.05
    alpha <- 0.36
    beta <- 0.21
    phi <- 0.9
    gamma <- 0.2
    l_0 <- 1
    b_0 <- 0.05
    s_ad <- c(-0.01, 0.01, 0.03, -0.03)
    s_m <- c(0.99, 1.01, 1.03, 0.97)
    h <- 1
    if (out) {
        epsilon <- 0.05 #Fraction of outliers
        K <- 20
        u_out <- unlist(lapply(rep(0, ts + 1), function(x)
                if (runif(1, 0, 1) < (1 - epsilon)) 0
                else rnorm(1, 0, K^2 * sigma^2)))
    } else {
        u_out <- rep(0, ts + 1)
    }
    e <- rnorm(ts + 1, 0, sigma^2)
    b <- rep(0, ts + 1)
    b[1] <- b_0
    l <- rep(0, ts + 1)
    l[1] <- l_0
    s <- rep(0, ts + m + 1)
    h_m <- floor((h - 1) %% m) + 1
    y <- rep(NA, ts + 1)
    u <- rep(0, ts + 1)
    y_pred <- rep(NA, ts + 1)
    first_char <- unlist(strsplit((type), split = ""))[1]

    if (type %in% c("ADN", "MDN", "ANN", "MNN", "AAN", "MAN")) {
        if (type %in% c("ANN", "MNN")) {
            #None Trend
            phi <- 0
            if (first_char  == "A") {
                for (t in seq(2, ts + 1)) {
                    u[[t]] <- l[[t - 1]]
                    l[[t]] <- l[[t - 1]] + alpha * e[[t]]
                }
                y <- u + e + u_out
            } else {
                for (t in seq(2, ts + 1)) {
                    u[[t]] <- l[[t - 1]]
                    l[[t]] <- l[[t - 1]] * (1 + alpha * e[[t]] * u[[t]])
                }
                y <- u * unlist(lapply((e + u_out), outlier_multi))
            }
        } else if (type %in% c("AAN", "MAN")) {
            # Additive Trend
            phi <- 1
            if (first_char  == "A") {
                for (t in seq(2, ts + 1)) {
                    u[[t]] <- l[[t - 1]] + b[[t - 1]]
                    l[[t]] <- l[[t - 1]] + b[[t - 1]] +  alpha * e[[t]]
                    b[[t]] <- b[[t - 1]] + alpha * beta * e[[t]]
                }
                y <- u + e + u_out
            } else {
                for (t in seq(2, ts + 1)) {
                    u[[t]] <- l[[t - 1]] + b[[t - 1]]
                    l[[t]] <- l[[t - 1]] + b[[t - 1]] +  alpha * e[[t]] * u[[t]]
                    b[[t]] <- b[[t - 1]] + alpha * beta * e[[t]] * u[[t]]
                }
                y <- u * unlist(lapply((e + u_out), outlier_multi))
            }
        } else {
            # Damped Trend
            if (first_char  == "A") {
                for (t in seq(2, ts + 1)) {
                    u[[t]] <- l[[t - 1]] + b[[t - 1]]
                    b[[t]] <- phi * b[[t - 1]] + alpha * beta * e[[t]]
                    l[[t]] <- l[[t - 1]] * b[[t - 1]] +  alpha * e[[t]]
                    }
                y <- u + e + u_out
            } else {
                for (t in seq(2, ts + 1)) {
                    u[[t]] <- l[[t - 1]] + b[[t - 1]]
                    b[[t]] <- phi * b[[t - 1]] + alpha * beta * e[[t]] * u[[t]]
                    l[[t]] <- l[[t - 1]] * b[[t - 1]] +  alpha * e[[t]] * u[[t]]
                }
                y <- u * unlist(lapply((e + u_out), outlier_multi))
            }
        }
        for (t in seq(2, ts + 1)) {
            l[[t]] <- alpha * y[[t]] +
             (1 - alpha) * (l[[t - 1]] + phi * b[[t - 1]])
            b[[t]] <- beta * (l[[t]] - l[[t - 1]]) +
             (1 - beta) * phi * b[[t - 1]]
            y_pred[t + h] <- l[[t]] + phi * b[[t]]
        }
    } else if (type %in% c("ANA", "MNA", "AAA", "MAA", "ADA", "MDA")) {
        s[1:4] <- s_ad
        #Underlying model
        if (type %in% c("ANA", "MNA")) {
            #None Trend
            phi <- 0
            if (first_char  == "A") { 
                for (t in seq(2, ts + 1)) {
                    t_offset <- t + m - 1
                    s[[t_offset]] <- s[[t_offset - m]] + gamma * e[[t]]
                    u[[t]] <- l[[t - 1]] + s[[t_offset - m]]
                    l[[t]] <- l[[t - 1]] + alpha * e[[t]]
                }
                y <- u + e + u_out
            } else {
                for (t in seq(2, ts + 1)) {
                        t_offset <- t + m - 1
                        u[[t]] <- l[[t - 1]] + s[[t_offset - m]]
                        s[[t_offset]] <- s[[t_offset - m]] +
                         gamma * e[[t]] * u[[t]]
                        l[[t]] <- l[[t - 1]] + alpha * e[[t]] * u[[t]]
                }
                y <- u * unlist(lapply((e + u_out), outlier_multi))
            }
        } else if (type %in% c("AAA", "MAA")) {
            # Additive Trend
            phi <- 1
            if (first_char  == "A") {
                for (t in seq(2, ts + 1)) {
                    t_offset <- t + m - 1
                    u[[t]] <- l[[t - 1]] + b[[t - 1]] + s[[t_offset - m]]
                    s[[t_offset]] <- s[[t_offset - m]] + gamma * e[[t]]
                    b[[t]] <- b[[t - 1]] + alpha * beta * e[[t]]
                    l[[t]] <- l[[t - 1]] + b[[t - 1]] + alpha * e[[t]]
                }
                y <- u + e + u_out
            } else {
                for (t in seq(2, ts + 1)) {
                    t_offset <- t + m - 1
                    u[[t]] <- l[[t - 1]] + b[[t - 1]] + s[[t_offset - m]]
                    s[[t_offset]] <- s[[t_offset - m]] + gamma * e[[t]] * u[[t]]
                    b[[t]] <- b[[t - 1]] + alpha * beta * e[[t]] * u[[t]]
                    l[[t]] <- l[[t - 1]] + b[[t - 1]] + alpha * e[[t]] * u[[t]]
                }
                y <- u * unlist(lapply((e + u_out), outlier_multi))
            }
        } else {
            #Damped Trend
            if (first_char  == "A") {
                for (t in seq(2, ts + 1)) {
                    t_offset <- t + m - 1
                    u[[t]] <- l[[t - 1]] + b[[t - 1]] + s[[t_offset - m]]
                    s[[t_offset]] <- s[[t_offset - m]] + gamma * e[[t]]
                    b[[t]] <- phi * b[[t - 1]] + alpha * beta * e[[t]]
                    l[[t]] <- l[[t - 1]] + b[[t - 1]] + alpha * e[[t]]
                }
                y <- u + e + u_out
            } else {
                for (t in seq(2, ts + 1)) {
                    t_offset <- t + m - 1
                    u[[t]] <- l[[t - 1]] + b[[t - 1]] + s[[t_offset - m]]
                    s[[t_offset]] <- s[[t_offset - m]] + gamma * e[[t]] * u[[t]]
                    b[[t]] <- phi * b[[t - 1]] + alpha * beta * e[[t]] * u[[t]]
                    l[[t]] <- l[[t - 1]] + b[[t - 1]] + alpha * e[[t]] * u[[t]]
                }
                y <- u * unlist(lapply((e + u_out), outlier_multi))
            }
        }
        #Simulating data
        for (t in seq(2, ts + 1)) {
            t_offset <- t + m - 1
            s[[t_offset]] <- gamma * (y[[t]] - l[[t - 1]] - phi * b[[t - 1]]) +
             (1 - gamma) * s[[t_offset - m]]
            l[[t]] <- alpha * (y[[t]] - s[[t_offset - m]]) +
             (1 - alpha) * (l[[t - 1]] + phi * b[[t - 1]])
            b[[t]] <- beta * (l[[t]] - l[[t - 1]]) +
             (1 - beta) * phi * b[[t - 1]]
            y_pred[t + h] <- l[[t]] + phi * b[[t]] + s[t_offset - m + h_m]
        }
    } else if (type %in% c("MNM", "MAM", "MDM")) {
        s[1:4] <- s_m
        if (type %in% c("MNM")) {
            phi <- 0
            for (t in seq(2, ts + 1)) {
                t_offset <- t + m - 1
                u[[t]] <- l[[t - 1]] * s[[t_offset - m]]
                s[[t_offset]] <- s[[t_offset - m]] +
                 (gamma * e[[t]] * u[[t]]) / l[[t - 1]]
                l[[t]] <- l[[t - 1]] +
                 (alpha * e[[t]] * u[[t]]) / s[[t_offset - m]]
            }
            y <- u * unlist(lapply((e + u_out), outlier_multi))
        } else if (type %in% c("MAM")) {
            phi <- 1
            for (t in seq(2, ts + 1)) {
                t_offset <- t + m - 1
                u[[t]] <- ((l[[t - 1]] + b[[t - 1]]) * s[[t_offset - m]])
                s[[t_offset]] <- s[[t_offset - m]] +
                 (gamma * e[[t]] * u[[t]]) / (l[[t - 1]] + b[[t - 1]])
                b[[t]] <- b[[t - 1]] +
                 (alpha * beta * e[[t]] * u[[t]]) / s[[t_offset - m]]
                l[[t]] <- l[[t - 1]] +
                 b[[t - 1]] + (alpha * e[[t]] * u[[t]]) / s[[t_offset - m]]
            }
            y <- u * unlist(lapply((e + u_out), outlier_multi))
        } else {
            for (t in seq(2, ts + 1)) {
                t_offset <- t + m - 1
                u[[t]] <- ((l[[t - 1]] + b[[t - 1]]) * s[[t_offset - m]])
                s[[t_offset]] <- s[[t_offset - m]] +
                 (gamma * e[[t]] * u[[t]]) / (l[[t - 1]] + b[[t - 1]])
                b[[t]] <- phi * b[[t - 1]] +
                 (alpha * beta * e[[t]] * u[[t]]) / s[[t_offset - m]]
                l[[t]] <- l[[t - 1]] + b[[t - 1]] +
                 (alpha * e[[t]] * u[[t]]) / s[[t_offset - m]]
            }
            y <- u * unlist(lapply((e + u_out), outlier_multi))
        }
        for (t in seq(1, ts)) {
            t_offset <- t + m - 1
            s[[t_offset]] <- gamma * (y[[t]] / (l[[t - 1]] +
             phi * b[[t - 1]])) + (1 - gamma) * s[[t_offset - m]]
            l[[t]] <- alpha * (y[[t]] / s[[t_offset - m]]) +
             (1 - alpha) * (l[[t - 1]] + phi * b[[t - 1]])
            b[[t]] <- beta * (l[[t]] - l[[t - 1]]) +
             (1 - beta) * phi * b[[t - 1]]
            y_pred[t + h] <- (l[[t]] + phi * b[[t]]) * s[t_offset - m + h_m]
        }
    } else {
        print("Error: Model is not initialized.")
    }
    return(list("y" = y, "y_pred" = y_pred, "model" = type))
}


get_best_psi <- function(x, model, scale.estimator = "mad", psifun = "Huber", damping = FALSE) {
    l <- c()
    a <- c("Huber", "Bisquare", "Hampel", "Welsh1", "Welsh2")
    for (i in a) {
        m_ahenao <- robets(x,
        model = model,
        scale.estimator = scale.estimator,
        psifun = i,
        damped = damping)
        l[[paste0(i)]] <- m_ahenao$robaicc
    }
    return(which.min(l))
}
simulate_results <- function(i = 1000, m = "ANN", damping = FALSE) {
    model_type <- m
    iterations <- i
    r <- vector("list", iterations)
    if (is.element(substr(m, 2, 2), c("D"))) {
        substr(model_type, 2, 2) <- "A"
        damping <- TRUE
    }
    for (n in seq(iterations)) {
        m.t <- model(model_type, out = FALSE)  # Outlier flag
        if (model_type %in% c("ADN", "MDN")) {
        time_series <- m.t$y[c(-1, -2)]
        } else {
        time_series <- m.t$y[c(-1)]
        }
        title <- get_best_psi(time_series,
         model = model_type,
         damping = damping)
        mf <- mean(time_series[1:11])
        ml <- mean(time_series[31:length(time_series)])
        mc <- abs(mf - ml) / mf
        r[[n]] <- c(n, names(title), model_type,
        sd(time_series), mf, ml, mc
        )
    }
    data <- do.call(rbind, r)
    colnames(data) <- c("iter", "loss", "model", "std", "mean_first",
     "mean_last", "mean_change")
    data <- transform(data, iter =  as.numeric(iter),
    std = as.numeric(std),
    mean_first = as.numeric(mean_first),
    mean_last = as.numeric(mean_last),
    mean_change = as.numeric(mean_change)
    )
    return(data)
}


# result <- data %>%
#             group_by(loss) %>%
#              summarise(std = median(std),
#               mean_first = median(mean_first),
#               mean_last = median(mean_last),
#               mean_change = median(mean_change),
#               model = n() / 10)
# "ANN", "MNN", "AAN", "MAN", "AAdN", "MAdN", "ANA", "MNA", "AAA", "MAA", "AAdA", "MAdA"

res <- list()
for (i in seq(10)) {
    data <- simulate_results(1000, "MNA")
    x <- data %>%
         group_by(loss) %>%
         summarise(std = median(std),
          mean_first = median(mean_first),
          mean_last = median(mean_last),
          mean_change = median(mean_change),
          model = n() / 10)
    res[[i]] <- x
}
binded <- bind_rows(res)
print(binded, n = 100)

## UNUSED UNDERLYING MODELS
# #AMN
# for (t in seq(2, ts + 1)) {
#     b[[t]] <- b[[t - 1]] + (alpha * beta * e[[t]]) / l[[t - 1]]
#     l[[t]] <- l[[t - 1]] * b[[t - 1]] +  alpha * e[[t]]
#     u[[t]] <- l[[t - 1]] * b[[t - 1]]
#     y[[t]] <- u[[t]] + e[[t]]
# }
# #AMA
# s[:3] <- s_ad
# for (t in seq(ts)) {
#     t_offset <- t + m - 1
#     s[[t_offset]] = s[[t_offset - m]] + gamma * e[[t]]
#     b[[t]] <- b[[t - 1]] + (alpha * beta * e[[t]]) / l[[t - 1]]
#     l[[t]] <- l[[t - 1]] * b[[t - 1]] + alpha * e[[t]]
#     y[[t]] <- l[[t - 1]] * b[[t - 1]] + s[[t_offset - m]] + e[[t]]
# }
# #ANM
# s[:3] <- s_m
# for (t in seq(ts)) {
#     t_offset <- t + m - 1
#     s[[t_offset]] = s[[t_offset - m]] + (gamma * e[[t]]) / l[[t - 1]]
#     l[[t]] <- l[[t - 1]] + (alpha * e[[t]]) / s[[t_offset - m]]
#     y[[t]] <- l[[t - 1]] * s[[t_offset - m]] + e[[t]]
# }
# #AAM
# s[:3] <- s_m
# for (t in seq(ts)) {
#     t_offset <- t + m - 1
#     s[[t_offset]] = s[[t_offset - m]] + (gamma * e[[t]]) / (l[[t - 1]] + b[[t - 1]])
#     b[[t]] <- b[[t - 1]] + (alpha * beta * e[[t]]) / s[[t_offset - m]]
#     l[[t]] <- l[[t - 1]] + b[[t - 1]] + (alpha * e[[t]]) / s[[t_offset - m]]
#     y[[t]] <- ((l[[t - 1]] + b[[t - 1]]) * s[[t_offset - m]]) + e[[t]]
# }
# #AMM
# s[:3] <- s_m
# for (t in seq(ts)) {
#     t_offset <- t + m - 1
#     s[[t_offset]] = s[[t_offset - m]] + (gamma * e[[t]]) / (l[[t - 1]] * b[[t - 1]])
#     b[[t]] <- b[[t - 1]] + (alpha * beta * e[[t]]) / (s[[t_offset - m]] * l[[t - 1]])
#     l[[t]] <- l[[t - 1]] * b[[t - 1]] + (alpha * e[[t]]) / s[[t_offset - m]]
#     y[[t]] <- ((l[[t - 1]] * b[[t - 1]]) * s[[t_offset - m]]) + e[[t]]
# }
# #AAdM
# s[:3] <- s_m
# for (t in seq(ts)) {
#     t_offset <- t + m - 1
#     s[[t_offset]] = s[[t_offset - m]] + (gamma * e[[t]]) / (l[[t - 1]] + b[[t - 1]])
#     b[[t]] <- phi * b[[t - 1]] + (alpha * beta * e[[t]]) / s[[t_offset - m]]
#     l[[t]] <- l[[t - 1]] + b[[t - 1]] + (alpha * e[[t]]) / s[[t_offset - m]]
#     y[[t]] <- ((l[[t - 1]] + b[[t - 1]]) * s[[t_offset - m]]) + e[[t]]
# }
