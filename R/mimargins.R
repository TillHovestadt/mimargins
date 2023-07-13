#' AME Calculation
#'
#' @param df The data frame which should be used for the regressions. Should be in long format and preferably not a mids-object. Whether mids-objects work has not yet been tested, however.
#' @param formula The Regression formula. Should preferably be a formula object. Must correspond to regression type specified in 'type'.
#' @param imp_var The variable in the data frame indicating the imputation (m). Defaults to ".imp".
#' @param ame Whether average marginal effects (AMEs) shall be calculated. Defaults to TRUE. If set to FALSE, a regression with the formula above is run for each imputation and pooled using Rubin's Rules. If set to true, AMEs are calculated for each imputation and then pooled using Rubin's rules.
#' @param variables Whether AMEs shall be calculated for a specific variable only. Currently only works for specifying one variable. Only has an effect if ame = T. Defaults to NULL. Then, AMEs are calculated.
#' @param at Whether AMEs shall be calculated at specified levels of another variable. Currently only works for specifying one variable. The number of levels may vary. Input has to take the form list(VARIABLENAME = c(level1, level2, level3 ...)). Only has an effect if ame = T. Defaults to NULL. Then, AMEs are averaged over all levels of all other variables.
#' @param linkfunction If type is set to 'glm', a linkfunction needs to be specified. Has no effect if type is set to other values.
#' @param type Specifies which type of regression to run. Currently supports 'lm' (Default), 'lmer' and 'glm'. Has to correspond to the specified formula.
#' @return Returns Dataframe of AMEs and SEs.
#' @export
mimargins <- function(df,
                      formula,
                      imp_var = ".imp",
                      type = "lm",
                      ame = T,
                      variables = NULL,
                      at = NULL,
                      linkfunction = NULL) {

  # for linear OLS regression ----
  if (type == "lm") {
    if (ame == T) {
      if (is.null(at)) {
        ifelse(!require(margins),
               install.packages("margins", dep = T),
               require(margins))
        require(margins)


        test  <- lm(formula = formula,
                    data = df[df[[imp_var]] == 1, ])
        a <- summary(margins::margins(test,
                                      variables = variables))

        imps <- length(unique(df[[imp_var]]))

        ame_coef <-
          data.frame(matrix(
            data = rep(x = NA, times = imps * length(a$factor)),
            ncol = length(a$factor)
          ))
        ame_se <-
          data.frame(matrix(
            data = rep(x = NA, times = imps * length(a$factor)),
            ncol = length(a$factor)
          ))

        colnames(ame_coef) <- a$factor
        colnames(ame_se) <- a$factor

        ame_coef$imputation <- c(1:imps)
        ame_se$imputation <- c(1:imps)

        for (imp in unique(df[[imp_var]])) {
          reg <- lm(formula = formula,
                    data = df[df[[imp_var]] == imp, ])

          regmargins <- summary(margins::margins(reg,
                                                 at = at,
                                                 variables = variables))

          for (eff in regmargins$factor) {
            ame_coef[[paste0(eff)]][ame_coef$imputation == imp] <-
              regmargins$AME[regmargins$factor == paste0(eff)]
            ame_se[[paste0(eff)]][ame_se$imputation == imp] <-
              regmargins$SE[regmargins$factor == paste0(eff)]

          }



        }

        # calculate between variance

        mean_ame <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        within_var <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        between_var <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        pooled_se <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))

        colnames(mean_ame) <- a$factor
        colnames(within_var) <- a$factor
        colnames(between_var) <- a$factor
        colnames(pooled_se) <- a$factor

        for (eff in a$factor) {
          mean_ame[[paste0(eff)]] <- mean(x = ame_coef[[paste0(eff)]])
          within_var[[paste0(eff)]] <-
            mean(x = ame_se[[paste0(eff)]] ^ 2)
          between_var[[paste0(eff)]] <-
            var(x = ame_coef[[paste0(eff)]])
          pooled_se[[paste0(eff)]] <-
            sqrt(within_var[[paste0(eff)]] + between_var[[paste0(eff)]] + (between_var[[paste0(eff)]] / imps))

        }

        pred <- data.frame(name = a$factor,
                           AME = NA,
                           SE = NA)

        for (eff in a$factor) {
          pred$AME[pred$name == paste0(eff)] <- mean_ame[[paste0(eff)]]
          pred$SE[pred$name == paste0(eff)] <-
            pooled_se[[paste0(eff)]]

        }

        return(pred)




      } else {
        ifelse(!require(margins),
               install.packages("margins", dep = T),
               require(margins))
        require(margins)


        test  <- lm(formula = formula,
                    data = df[df[[imp_var]] == 1, ])
        a <- summary(margins::margins(test,
                                      variables = variables,
                                      at = at))

        imps <- length(unique(df[[imp_var]]))

        ame_coef <-
          data.frame(matrix(
            data = rep(x = NA, times = imps * length(a$factor)),
            ncol = length(a$factor)
          ))
        ame_se <-
          data.frame(matrix(
            data = rep(x = NA, times = imps * length(a$factor)),
            ncol = length(a$factor)
          ))

        colnames(ame_coef) <- paste0(a$factor,
                                     "_at_",
                                     names(x = at),
                                     "_=_",
                                     a[[names(x = at)]])
        colnames(ame_se) <- paste0(a$factor,
                                   "_at_",
                                   names(x = at),
                                   "_=_",
                                   a[[names(x = at)]])

        ame_coef$imputation <- c(1:imps)
        ame_se$imputation <- c(1:imps)

        for (imp in unique(df[[imp_var]])) {
          reg <- lm(formula = formula,
                    data = df[df[[imp_var]] == imp, ])

          regmargins <- summary(margins::margins(reg,
                                                 at = at,
                                                 variables = variables))

          columnnames <- paste0(regmargins$factor,
                                "_at_",
                                names(x = at),
                                "_=_",
                                regmargins[[names(x = at)]])



          for (eff in unique(regmargins$factor)) {
            for (value in regmargins[[names(at)]]) {
              ame_coef[[paste0(eff,
                               "_at_",
                               names(x = at),
                               "_=_",
                               value)]][ame_coef$imputation == imp] <-
                regmargins$AME[(regmargins$factor == paste0(eff) &
                                  regmargins[[names(at)]] == value)]

              ame_se[[paste0(eff,
                             "_at_",
                             names(x = at),
                             "_=_",
                             value)]][ame_se$imputation == imp] <-
                regmargins$SE[(regmargins$factor == paste0(eff) &
                                 regmargins[[names(at)]] == value)]




            }


          }



        }

        # calculate between variance

        mean_ame <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        within_var <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        between_var <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        pooled_se <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))

        colnames(mean_ame) <- paste0(a$factor,
                                     "_at_",
                                     names(x = at),
                                     "_=_",
                                     a[[names(x = at)]])
        colnames(within_var) <- paste0(a$factor,
                                       "_at_",
                                       names(x = at),
                                       "_=_",
                                       a[[names(x = at)]])
        colnames(between_var) <- paste0(a$factor,
                                        "_at_",
                                        names(x = at),
                                        "_=_",
                                        a[[names(x = at)]])
        colnames(pooled_se) <- paste0(a$factor,
                                      "_at_",
                                      names(x = at),
                                      "_=_",
                                      a[[names(x = at)]])

        for (eff in colnames(mean_ame)) {
          mean_ame[[paste0(eff)]] <- mean(x = ame_coef[[paste0(eff)]])
          within_var[[paste0(eff)]] <-
            mean(x = ame_se[[paste0(eff)]] ^ 2)
          between_var[[paste0(eff)]] <-
            var(x = ame_coef[[paste0(eff)]])
          pooled_se[[paste0(eff)]] <-
            sqrt(within_var[[paste0(eff)]] + between_var[[paste0(eff)]] + (between_var[[paste0(eff)]] / imps))

        }

        pred <- data.frame(
          name = paste0(a$factor,
                        "_at_",
                        names(x = at),
                        "_=_",
                        a[[names(x = at)]]),
          AME = NA,
          SE = NA,
          at_level = a[[names(x = at)]]
        )

        for (eff in pred$name) {
          pred$AME[pred$name == paste0(eff)] <- mean_ame[[paste0(eff)]]
          pred$SE[pred$name == paste0(eff)] <-
            pooled_se[[paste0(eff)]]

        }

        return(pred)






      }

    }

    if (ame == F) {
      test  <- lm(formula = formula,
                  data = df[df[[imp_var]] == 1, ])
      a <- summary(test)
      a <- as.data.frame(a$coefficients)

      varnames <- rownames(a)

      imps <- length(unique(df[[imp_var]]))

      res_coef <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(varnames)),
          ncol = length(varnames)
        ))
      res_se <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(varnames)),
          ncol = length(varnames)
        ))

      colnames(res_coef) <- varnames
      colnames(res_se) <- varnames

      res_coef$imputation <- c(1:imps)
      res_se$imputation <- c(1:imps)

      res_coef$R2 <- NA

      for (imp in unique(df[[imp_var]])) {
        reg <- summary(lm(formula = formula,
                          data = df[df[[imp_var]] == imp, ]))

        for (eff in rownames(reg$coefficients)) {
          res_coef[[paste0(eff)]][res_coef$imputation == imp] <-
            reg$coefficients[rownames(reg$coefficients) == paste0(eff)][1]
          res_se[[paste0(eff)]][res_se$imputation == imp] <-
            reg$coefficients[rownames(reg$coefficients) == paste0(eff)][2]
          res_coef[["R2"]][res_coef$imputation == imp] <-
            reg$adj.r.squared
        }



      }

      # calculate between variance

      mean_coef <-
        data.frame(matrix(
          data = rep(x = NA, times = length(varnames)),
          ncol = length(varnames)
        ))
      within_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(varnames)),
          ncol = length(varnames)
        ))
      between_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(varnames)),
          ncol = length(varnames)
        ))
      pooled_se <-
        data.frame(matrix(
          data = rep(x = NA, times = length(varnames)),
          ncol = length(varnames)
        ))

      mean_r2 <- mean(res_coef$R2)

      colnames(mean_coef) <- varnames
      colnames(within_var) <- varnames
      colnames(between_var) <- varnames
      colnames(pooled_se) <- varnames

      for (eff in varnames) {
        mean_coef[[paste0(eff)]] <- mean(x = res_coef[[paste0(eff)]])
        within_var[[paste0(eff)]] <-
          mean(x = res_se[[paste0(eff)]] ^ 2)
        between_var[[paste0(eff)]] <-
          var(x = res_coef[[paste0(eff)]])
        pooled_se[[paste0(eff)]] <-
          sqrt(within_var[[paste0(eff)]] + between_var[[paste0(eff)]] + (between_var[[paste0(eff)]] / imps))


      }

      res <- data.frame(
        name = varnames,
        AME = NA,
        SE = NA,
        R2 = mean_r2
      )

      for (eff in varnames) {
        res$AME[res$name == paste0(eff)] <- mean_coef[[paste0(eff)]]
        res$SE[res$name == paste0(eff)] <- pooled_se[[paste0(eff)]]

      }

      return(res)




    }

  }

  ## for hierarchical linear regressions ----

  if (type == "lmer") {
    ifelse(!require(lme4),
           install.packages("lme4", dep = T),
           require(lme4))
    require(lme4)


    if (ame == T) {
      if (is.null(at)) {
        ifelse(!require(margins),
               install.packages("margins", dep = T),
               require(margins))
        require(margins)


        test  <- lmer(formula = formula,
                      data = df[df[[imp_var]] == 1, ])
        a <- summary(margins::margins(test,
                                      variables = variables))

        imps <- length(unique(df[[imp_var]]))

        ame_coef <-
          data.frame(matrix(
            data = rep(x = NA, times = imps * length(a$factor)),
            ncol = length(a$factor)
          ))
        ame_se <-
          data.frame(matrix(
            data = rep(x = NA, times = imps * length(a$factor)),
            ncol = length(a$factor)
          ))

        colnames(ame_coef) <- a$factor
        colnames(ame_se) <- a$factor

        ame_coef$imputation <- c(1:imps)
        ame_se$imputation <- c(1:imps)

        for (imp in unique(df[[imp_var]])) {
          reg <- lmer(formula = formula,
                      data = df[df[[imp_var]] == imp, ])

          regmargins <- summary(margins::margins(reg,
                                                 at = at,
                                                 variables = variables))

          for (eff in regmargins$factor) {
            ame_coef[[paste0(eff)]][ame_coef$imputation == imp] <-
              regmargins$AME[regmargins$factor == paste0(eff)]
            ame_se[[paste0(eff)]][ame_se$imputation == imp] <-
              regmargins$SE[regmargins$factor == paste0(eff)]

          }



        }

        # calculate between variance

        mean_ame <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        within_var <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        between_var <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        pooled_se <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))

        colnames(mean_ame) <- a$factor
        colnames(within_var) <- a$factor
        colnames(between_var) <- a$factor
        colnames(pooled_se) <- a$factor

        for (eff in a$factor) {
          mean_ame[[paste0(eff)]] <- mean(x = ame_coef[[paste0(eff)]])
          within_var[[paste0(eff)]] <-
            mean(x = ame_se[[paste0(eff)]] ^ 2)
          between_var[[paste0(eff)]] <-
            var(x = ame_coef[[paste0(eff)]])
          pooled_se[[paste0(eff)]] <-
            sqrt(within_var[[paste0(eff)]] + between_var[[paste0(eff)]] + (between_var[[paste0(eff)]] / imps))

        }

        pred <- data.frame(name = a$factor,
                           AME = NA,
                           SE = NA)

        for (eff in a$factor) {
          pred$AME[pred$name == paste0(eff)] <- mean_ame[[paste0(eff)]]
          pred$SE[pred$name == paste0(eff)] <-
            pooled_se[[paste0(eff)]]

        }

        return(pred)




      } else {
        ifelse(!require(margins),
               install.packages("margins", dep = T),
               require(margins))
        require(margins)


        test  <- lmer(formula = formula,
                      data = df[df[[imp_var]] == 1, ])
        a <- summary(margins::margins(test,
                                      variables = variables,
                                      at = at))

        imps <- length(unique(df[[imp_var]]))

        ame_coef <-
          data.frame(matrix(
            data = rep(x = NA, times = imps * length(a$factor)),
            ncol = length(a$factor)
          ))
        ame_se <-
          data.frame(matrix(
            data = rep(x = NA, times = imps * length(a$factor)),
            ncol = length(a$factor)
          ))

        colnames(ame_coef) <- paste0(a$factor,
                                     "_at_",
                                     names(x = at),
                                     "_=_",
                                     a[[names(x = at)]])
        colnames(ame_se) <- paste0(a$factor,
                                   "_at_",
                                   names(x = at),
                                   "_=_",
                                   a[[names(x = at)]])

        ame_coef$imputation <- c(1:imps)
        ame_se$imputation <- c(1:imps)

        for (imp in unique(df[[imp_var]])) {
          reg <- lmer(formula = formula,
                      data = df[df[[imp_var]] == imp, ])

          regmargins <- summary(margins::margins(reg,
                                                 at = at,
                                                 variables = variables))

          columnnames <- paste0(regmargins$factor,
                                "_at_",
                                names(x = at),
                                "_=_",
                                regmargins[[names(x = at)]])



          for (eff in unique(regmargins$factor)) {
            for (value in regmargins[[names(at)]]) {
              ame_coef[[paste0(eff,
                               "_at_",
                               names(x = at),
                               "_=_",
                               value)]][ame_coef$imputation == imp] <-
                regmargins$AME[(regmargins$factor == paste0(eff) &
                                  regmargins[[names(at)]] == value)]

              ame_se[[paste0(eff,
                             "_at_",
                             names(x = at),
                             "_=_",
                             value)]][ame_se$imputation == imp] <-
                regmargins$SE[(regmargins$factor == paste0(eff) &
                                 regmargins[[names(at)]] == value)]




            }


          }



        }

        # calculate between variance

        mean_ame <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        within_var <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        between_var <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        pooled_se <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))

        colnames(mean_ame) <- paste0(a$factor,
                                     "_at_",
                                     names(x = at),
                                     "_=_",
                                     a[[names(x = at)]])
        colnames(within_var) <- paste0(a$factor,
                                       "_at_",
                                       names(x = at),
                                       "_=_",
                                       a[[names(x = at)]])
        colnames(between_var) <- paste0(a$factor,
                                        "_at_",
                                        names(x = at),
                                        "_=_",
                                        a[[names(x = at)]])
        colnames(pooled_se) <- paste0(a$factor,
                                      "_at_",
                                      names(x = at),
                                      "_=_",
                                      a[[names(x = at)]])

        for (eff in colnames(mean_ame)) {
          mean_ame[[paste0(eff)]] <- mean(x = ame_coef[[paste0(eff)]])
          within_var[[paste0(eff)]] <-
            mean(x = ame_se[[paste0(eff)]] ^ 2)
          between_var[[paste0(eff)]] <-
            var(x = ame_coef[[paste0(eff)]])
          pooled_se[[paste0(eff)]] <-
            sqrt(within_var[[paste0(eff)]] + between_var[[paste0(eff)]] + (between_var[[paste0(eff)]] / imps))

        }

        pred <- data.frame(
          name = paste0(a$factor,
                        "_at_",
                        names(x = at),
                        "_=_",
                        a[[names(x = at)]]),
          AME = NA,
          SE = NA,
          at_level = a[[names(x = at)]]
        )

        for (eff in pred$name) {
          pred$AME[pred$name == paste0(eff)] <- mean_ame[[paste0(eff)]]
          pred$SE[pred$name == paste0(eff)] <-
            pooled_se[[paste0(eff)]]

        }

        return(pred)






      }

    }

    if (ame == F) {
      test  <- lmer(formula = formula,
                    data = df[df[[imp_var]] == 1, ])
      a <- summary(test)
      a <- as.data.frame(a$coefficients)

      varnames <- rownames(a)

      imps <- length(unique(df[[imp_var]]))

      res_coef <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(varnames)),
          ncol = length(varnames)
        ))
      res_se <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(varnames)),
          ncol = length(varnames)
        ))

      colnames(res_coef) <- varnames
      colnames(res_se) <- varnames

      res_coef$imputation <- c(1:imps)
      res_se$imputation <- c(1:imps)

      res_coef$R2 <- NA

      for (imp in unique(df[[imp_var]])) {
        reg <- summary(lmer(formula = formula,
                            data = df[df[[imp_var]] == imp, ]))

        for (eff in rownames(reg$coefficients)) {
          res_coef[[paste0(eff)]][res_coef$imputation == imp] <-
            reg$coefficients[rownames(reg$coefficients) == paste0(eff)][1]
          res_se[[paste0(eff)]][res_se$imputation == imp] <-
            reg$coefficients[rownames(reg$coefficients) == paste0(eff)][2]
          res_coef[["R2"]][res_coef$imputation == imp] <-
            reg$adj.r.squared
        }



      }

      # calculate between variance

      mean_coef <-
        data.frame(matrix(
          data = rep(x = NA, times = length(varnames)),
          ncol = length(varnames)
        ))
      within_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(varnames)),
          ncol = length(varnames)
        ))
      between_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(varnames)),
          ncol = length(varnames)
        ))
      pooled_se <-
        data.frame(matrix(
          data = rep(x = NA, times = length(varnames)),
          ncol = length(varnames)
        ))

      mean_r2 <- mean(res_coef$R2)

      colnames(mean_coef) <- varnames
      colnames(within_var) <- varnames
      colnames(between_var) <- varnames
      colnames(pooled_se) <- varnames

      for (eff in varnames) {
        mean_coef[[paste0(eff)]] <- mean(x = res_coef[[paste0(eff)]])
        within_var[[paste0(eff)]] <-
          mean(x = res_se[[paste0(eff)]] ^ 2)
        between_var[[paste0(eff)]] <-
          var(x = res_coef[[paste0(eff)]])
        pooled_se[[paste0(eff)]] <-
          sqrt(within_var[[paste0(eff)]] + between_var[[paste0(eff)]] + (between_var[[paste0(eff)]] / imps))


      }

      res <- data.frame(
        name = varnames,
        AME = NA,
        SE = NA,
        R2 = mean_r2
      )

      for (eff in varnames) {
        res$AME[res$name == paste0(eff)] <- mean_coef[[paste0(eff)]]
        res$SE[res$name == paste0(eff)] <- pooled_se[[paste0(eff)]]

      }

      return(res)




    }

  }

  # for glm models ----

  if (type == "glm") {
    stopifnot(!is.null(linkfunction))



    if (ame == T) {
      if (is.null(at)) {
        ifelse(!require(margins),
               install.packages("margins", dep = T),
               require(margins))
        require(margins)


        test  <- glm(formula = formula,
                     data = df[df[[imp_var]] == 1, ],
                     family = linkfunction)
        a <- summary(margins::margins(test,
                                      variables = variables))

        imps <- length(unique(df[[imp_var]]))

        ame_coef <-
          data.frame(matrix(
            data = rep(x = NA, times = imps * length(a$factor)),
            ncol = length(a$factor)
          ))
        ame_se <-
          data.frame(matrix(
            data = rep(x = NA, times = imps * length(a$factor)),
            ncol = length(a$factor)
          ))

        colnames(ame_coef) <- a$factor
        colnames(ame_se) <- a$factor

        ame_coef$imputation <- c(1:imps)
        ame_se$imputation <- c(1:imps)

        for (imp in unique(df[[imp_var]])) {
          reg <- glm(formula = formula,
                     data = df[df[[imp_var]] == imp, ],
                     family = linkfunction)

          regmargins <- summary(margins::margins(reg,
                                                 at = at,
                                                 variables = variables))

          for (eff in regmargins$factor) {
            ame_coef[[paste0(eff)]][ame_coef$imputation == imp] <-
              regmargins$AME[regmargins$factor == paste0(eff)]
            ame_se[[paste0(eff)]][ame_se$imputation == imp] <-
              regmargins$SE[regmargins$factor == paste0(eff)]

          }



        }

        # calculate between variance

        mean_ame <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        within_var <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        between_var <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        pooled_se <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))

        colnames(mean_ame) <- a$factor
        colnames(within_var) <- a$factor
        colnames(between_var) <- a$factor
        colnames(pooled_se) <- a$factor

        for (eff in a$factor) {
          mean_ame[[paste0(eff)]] <- mean(x = ame_coef[[paste0(eff)]])
          within_var[[paste0(eff)]] <-
            mean(x = ame_se[[paste0(eff)]] ^ 2)
          between_var[[paste0(eff)]] <-
            var(x = ame_coef[[paste0(eff)]])
          pooled_se[[paste0(eff)]] <-
            sqrt(within_var[[paste0(eff)]] + between_var[[paste0(eff)]] + (between_var[[paste0(eff)]] / imps))

        }

        pred <- data.frame(name = a$factor,
                           AME = NA,
                           SE = NA)

        for (eff in a$factor) {
          pred$AME[pred$name == paste0(eff)] <- mean_ame[[paste0(eff)]]
          pred$SE[pred$name == paste0(eff)] <-
            pooled_se[[paste0(eff)]]

        }

        return(pred)




      } else {
        ifelse(!require(margins),
               install.packages("margins", dep = T),
               require(margins))
        require(margins)


        test  <- glm(formula = formula,
                     data = df[df[[imp_var]] == 1, ],
                     family = linkfunction)
        a <- summary(margins::margins(test,
                                      variables = variables,
                                      at = at))

        imps <- length(unique(df[[imp_var]]))

        ame_coef <-
          data.frame(matrix(
            data = rep(x = NA, times = imps * length(a$factor)),
            ncol = length(a$factor)
          ))
        ame_se <-
          data.frame(matrix(
            data = rep(x = NA, times = imps * length(a$factor)),
            ncol = length(a$factor)
          ))

        colnames(ame_coef) <- paste0(a$factor,
                                     "_at_",
                                     names(x = at),
                                     "_=_",
                                     a[[names(x = at)]])
        colnames(ame_se) <- paste0(a$factor,
                                   "_at_",
                                   names(x = at),
                                   "_=_",
                                   a[[names(x = at)]])

        ame_coef$imputation <- c(1:imps)
        ame_se$imputation <- c(1:imps)

        for (imp in unique(df[[imp_var]])) {
          reg <- glm(formula = formula,
                     data = df[df[[imp_var]] == imp, ],
                     family = linkfunction)

          regmargins <- summary(margins::margins(reg,
                                                 at = at,
                                                 variables = variables))

          columnnames <- paste0(regmargins$factor,
                                "_at_",
                                names(x = at),
                                "_=_",
                                regmargins[[names(x = at)]])



          for (eff in unique(regmargins$factor)) {
            for (value in regmargins[[names(at)]]) {
              ame_coef[[paste0(eff,
                               "_at_",
                               names(x = at),
                               "_=_",
                               value)]][ame_coef$imputation == imp] <-
                regmargins$AME[(regmargins$factor == paste0(eff) &
                                  regmargins[[names(at)]] == value)]

              ame_se[[paste0(eff,
                             "_at_",
                             names(x = at),
                             "_=_",
                             value)]][ame_se$imputation == imp] <-
                regmargins$SE[(regmargins$factor == paste0(eff) &
                                 regmargins[[names(at)]] == value)]




            }


          }



        }

        # calculate between variance

        mean_ame <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        within_var <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        between_var <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))
        pooled_se <-
          data.frame(matrix(
            data = rep(x = NA, times = length(a$factor)),
            ncol = length(a$factor)
          ))

        colnames(mean_ame) <- paste0(a$factor,
                                     "_at_",
                                     names(x = at),
                                     "_=_",
                                     a[[names(x = at)]])
        colnames(within_var) <- paste0(a$factor,
                                       "_at_",
                                       names(x = at),
                                       "_=_",
                                       a[[names(x = at)]])
        colnames(between_var) <- paste0(a$factor,
                                        "_at_",
                                        names(x = at),
                                        "_=_",
                                        a[[names(x = at)]])
        colnames(pooled_se) <- paste0(a$factor,
                                      "_at_",
                                      names(x = at),
                                      "_=_",
                                      a[[names(x = at)]])

        for (eff in colnames(mean_ame)) {
          mean_ame[[paste0(eff)]] <- mean(x = ame_coef[[paste0(eff)]])
          within_var[[paste0(eff)]] <-
            mean(x = ame_se[[paste0(eff)]] ^ 2)
          between_var[[paste0(eff)]] <-
            var(x = ame_coef[[paste0(eff)]])
          pooled_se[[paste0(eff)]] <-
            sqrt(within_var[[paste0(eff)]] + between_var[[paste0(eff)]] + (between_var[[paste0(eff)]] / imps))

        }

        pred <- data.frame(
          name = paste0(a$factor,
                        "_at_",
                        names(x = at),
                        "_=_",
                        a[[names(x = at)]]),
          AME = NA,
          SE = NA,
          at_level = a[[names(x = at)]]
        )

        for (eff in pred$name) {
          pred$AME[pred$name == paste0(eff)] <- mean_ame[[paste0(eff)]]
          pred$SE[pred$name == paste0(eff)] <-
            pooled_se[[paste0(eff)]]

        }

        return(pred)






      }

    }

    if (ame == F) {
      test  <- glm(formula = formula,
                   data = df[df[[imp_var]] == 1, ],
                   family = linkfunction)
      a <- summary(test)
      a <- as.data.frame(a$coefficients)

      varnames <- rownames(a)

      imps <- length(unique(df[[imp_var]]))

      res_coef <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(varnames)),
          ncol = length(varnames)
        ))
      res_se <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(varnames)),
          ncol = length(varnames)
        ))

      colnames(res_coef) <- varnames
      colnames(res_se) <- varnames

      res_coef$imputation <- c(1:imps)
      res_se$imputation <- c(1:imps)

      for (imp in unique(df[[imp_var]])) {
        reg <- summary(glm(
          formula = formula,
          data = df[df[[imp_var]] == imp, ],
          family = linkfunction
        ))

        for (eff in rownames(reg$coefficients)) {
          res_coef[[paste0(eff)]][res_coef$imputation == imp] <-
            reg$coefficients[rownames(reg$coefficients) == paste0(eff)][1]
          res_se[[paste0(eff)]][res_se$imputation == imp] <-
            reg$coefficients[rownames(reg$coefficients) == paste0(eff)][2]
        }



      }

      # calculate between variance

      mean_coef <-
        data.frame(matrix(
          data = rep(x = NA, times = length(varnames)),
          ncol = length(varnames)
        ))
      within_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(varnames)),
          ncol = length(varnames)
        ))
      between_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(varnames)),
          ncol = length(varnames)
        ))
      pooled_se <-
        data.frame(matrix(
          data = rep(x = NA, times = length(varnames)),
          ncol = length(varnames)
        ))

      colnames(mean_coef) <- varnames
      colnames(within_var) <- varnames
      colnames(between_var) <- varnames
      colnames(pooled_se) <- varnames

      for (eff in varnames) {
        mean_coef[[paste0(eff)]] <- mean(x = res_coef[[paste0(eff)]])
        within_var[[paste0(eff)]] <-
          mean(x = res_se[[paste0(eff)]] ^ 2)
        between_var[[paste0(eff)]] <-
          var(x = res_coef[[paste0(eff)]])
        pooled_se[[paste0(eff)]] <-
          sqrt(within_var[[paste0(eff)]] + between_var[[paste0(eff)]] + (between_var[[paste0(eff)]] / imps))


      }

      res <- data.frame(
        name = varnames,
        AME = NA,
        SE = NA
      )

      for (eff in varnames) {
        res$AME[res$name == paste0(eff)] <- mean_coef[[paste0(eff)]]
        res$SE[res$name == paste0(eff)] <- pooled_se[[paste0(eff)]]

      }

      return(res)




    }



  }




}
