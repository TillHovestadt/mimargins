#' Predicted Values Calculation
#'
#' @param df The data frame which should be used for the regressions. Should be in long format and preferably not a mids-object. Whether mids-objects work has not yet been tested, however.
#' @param formula The Regression formula. Should preferably be a formula object. Must correspond to regression type specified in 'type'.
#' @param imp_var The variable in the data frame indicating the imputation (m). Defaults to ".imp".
#' @param id_var The variable in the data frame differentiating the groups.
#' @param at Whether Predicted Values shall be calculated at specified levels of another variable. Currently only works for specifying one variable. The number of levels may vary. Input has to take the form list(VARIABLENAME = c(level1, level2, level3 ...)). Defaults to NULL. Then, Predicted Values are averaged over all levels of all other variables.
#' @param linkfunction If type is set to 'glm', a linkfunction needs to be specified. Has no effect if type is set to other values.
#' @param type Specifies which type of regression to run. Currently supports 'lm' (Default), 'lmer' and 'glm'. Has to correspond to the specified formula.
#' @param variables Specifies the variables of interest for which pooled estimates shall be computed. Only applies for type == "lmer". In all other cases, all variables in the regression model are taken as a default.
#' @param bootSeed For type = "lmer" only. Indicated whether a specific seed for the bootstraping of SEs shall be drawn. Otherwise, random seeds are drawn.
#' @param bootSims Number of Bootstrap Simulations to utilize. Defaults to 1000.
#' @return Returns Dataframe of Predicted Values and SEs.
#' @export
mipredictions <- function(df,
                          formula,
                          imp_var = ".imp",
                          id_var,
                          type = "lm",
                          at = NULL,
                          linkfunction = NULL,
                          variables = NULL,
                          bootSeed = NULL,
                          bootSims = 1000) {
  # for linear OLS regression ----
  if (type == "lm") {
    if (is.null(at)) {
      ifelse(!require(prediction),
             install.packages("prediction", dep = T),
             require(prediction))
      require(prediction)

      test  <- lm(formula = formula,
                  data = df[df[[imp_var]] == 1, ])

      a <- prediction::prediction(model = test,
                                  at = at,
                                  calculate_se = T)

      N <- length(a$fitted)
      imps <- length(unique(df[[imp_var]]))
      ids <- unique(df[[paste0(id_var)]])

      for (col in colnames(test$model)) {
        assign(x = paste0(col),
               value = data.frame(matrix(
                 data = rep(x = NA, times = imps * length(a$fitted)),
                 ncol = length(a$fitted)
               )))

        d <- get(paste0(col))
        colnames(d) <- a[[paste0(id_var)]]
        d$imputation <- c(1:imps)

        assign(x = paste0(col),
               value = d)

      }

      pred_coef <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(a$fitted)),
          ncol = length(a$fitted)
        ))
      pred_se <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(a$se.fitted)),
          ncol = length(a$se.fitted)
        ))

      colnames(pred_coef) <- a[[paste0(id_var)]]
      colnames(pred_se) <- a[[paste0(id_var)]]

      pred_coef$imputation <- c(1:imps)
      pred_se$imputation <- c(1:imps)

      for (imp in unique(df[[imp_var]])) {
        reg <- lm(formula = formula,
                  data = df[df[[imp_var]] == imp, ])

        regpred <- prediction::prediction(model = reg,
                                          at = at,
                                          calculate_se = T)

        regmodel <- reg$model
        regmodel$idvariable <- regpred[[paste0(id_var)]]

        for (id in ids) {
          pred_coef[[paste0(id)]][pred_coef$imputation == imp] <-
            regpred$fitted[regpred[[paste0(id_var)]] == id]
          pred_se[[paste0(id)]][pred_se$imputation == imp] <-
            regpred$se.fitted[regpred[[paste0(id_var)]] == id]

          for (col in colnames(reg$model)) {
            d <- get(paste0(col))

            d[[paste0(id)]][get(paste0(col))[["imputation"]] == imp] <-
              regmodel[[paste0(col)]][regmodel$idvariable == id]

            assign(x = paste0(col),
                   value = d)

          }

        }



      }

      # calculate between variance

      mean_pred <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      within_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      between_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      pooled_se <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))

      colnames(mean_pred) <- a[[paste0(id_var)]]
      colnames(within_var) <- a[[paste0(id_var)]]
      colnames(between_var) <- a[[paste0(id_var)]]
      colnames(pooled_se) <- a[[paste0(id_var)]]

      for (col in colnames(test$model)) {
        assign(x = paste0("mean_",
                          col),
               value = data.frame(matrix(
                 data = rep(x = NA, times = length(a$fitted)),
                 ncol = length(a$fitted)
               )))

        d <- get(paste0("mean_",
                        col))
        colnames(d) <- a[[paste0(id_var)]]

        assign(x = paste0("mean_",
                          col),
               value = d)

      }

      for (id in ids) {
        mean_pred[[paste0(id)]] <- mean(x = pred_coef[[paste0(id)]])
        within_var[[paste0(id)]] <-
          mean(x = pred_se[[paste0(id)]] ^ 2)
        between_var[[paste0(id)]] <-
          var(x = pred_coef[[paste0(id)]])
        pooled_se[[paste0(id)]] <-
          sqrt(within_var[[paste0(id)]] + between_var[[paste0(id)]] + (between_var[[paste0(id)]] / imps))

        for (col in colnames(test$model)) {
          d <- get(paste0("mean_",
                          col))

          d[[paste0(id)]] <-
            mean(x = get(paste0(col))[[paste0(id)]])

          assign(x = paste0("mean_",
                            col),
                 value = d)


        }

      }

      pred <- data.frame(
        pred = rep(x = NA,
                   times = length(ids)),
        pred_se = rep(x = NA,
                      times = length(ids))
      )
      pred[[paste0(id_var)]] <- ids

      for (col in colnames(test$model)) {
        pred[[paste0(col)]] <- NA



      }

      for (id in ids) {
        pred$pred[pred[[paste0(id_var)]] == id] <- mean_pred[[paste0(id)]]
        pred$pred_se[pred[[paste0(id_var)]] == id] <-
          pooled_se[[paste0(id)]]

        for (col in colnames(test$model)) {
          pred[[paste0(col)]][pred[[paste0(id_var)]] == id] <-
            get(paste0("mean_",
                       col))[[paste0(id)]]

        }

      }

      return(pred)




    } else {
      # for specified at() levels


      ifelse(!require(prediction),
             install.packages("prediction", dep = T),
             require(prediction))
      require(prediction)

      test  <- lm(formula = formula,
                  data = df[df[[imp_var]] == 1, ])

      a <- prediction::prediction(model = test,
                                  at = at,
                                  calculate_se = T)

      N <- length(a$fitted)
      imps <- length(unique(df[[imp_var]]))
      ids <- unique(df[[paste0(id_var)]])

      for (col in colnames(test$model)) {
        assign(x = paste0(col),
               value = data.frame(matrix(
                 data = rep(x = NA, times = imps * length(a$fitted)),
                 ncol = length(a$fitted)
               )))

        d <- get(paste0(col))
        colnames(d) <- paste0(a[[paste0(id_var)]],
                              "_at_",
                              names(x = at),
                              "_=_",
                              a[[names(x = at)]])
        d$imputation <- c(1:imps)

        assign(x = paste0(col),
               value = d)

      }

      pred_coef <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(a$fitted)),
          ncol = length(a$fitted)
        ))
      pred_se <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(a$se.fitted)),
          ncol = length(a$se.fitted)
        ))

      colnames(pred_coef) <- paste0(a[[paste0(id_var)]],
                                    "_at_",
                                    names(x = at),
                                    "_=_",
                                    a[[names(x = at)]])
      colnames(pred_se) <- paste0(a[[paste0(id_var)]],
                                  "_at_",
                                  names(x = at),
                                  "_=_",
                                  a[[names(x = at)]])

      pred_coef$imputation <- c(1:imps)
      pred_se$imputation <- c(1:imps)

      for (imp in unique(df[[imp_var]])) {
        reg <- lm(formula = formula,
                  data = df[df[[imp_var]] == imp, ])

        regpred <- prediction::prediction(model = reg,
                                          at = at,
                                          calculate_se = T)

        regmodel <- reg$model
        regmodel$idvariable <- unique(regpred[[paste0(id_var)]])

        for (id in ids) {
          for (value in regpred[[names(at)]]) {
            pred_coef[[paste0(id,
                              "_at_",
                              names(x = at),
                              "_=_",
                              value)]][pred_coef$imputation == imp] <-
              regpred$fitted[(regpred[[paste0(id_var)]] == id &
                                regpred[[names(at)]] == value)]

            pred_se[[paste0(id,
                            "_at_",
                            names(x = at),
                            "_=_",
                            value)]][pred_se$imputation == imp] <-
              regpred$se.fitted[(regpred[[paste0(id_var)]] == id &
                                   regpred[[names(at)]] == value)]

            for (col in colnames(reg$model)) {
              d <- get(paste0(col))

              d[[paste0(id,
                        "_at_",
                        names(x = at),
                        "_=_",
                        value)]][get(paste0(col))[["imputation"]] == imp] <-
                regmodel[[paste0(col)]][regmodel$idvariable == id]

              assign(x = paste0(col),
                     value = d)

            }


          }


        }



      }

      # calculate between variance

      mean_pred <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      within_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      between_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      pooled_se <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))

      colnames(mean_pred) <- paste0(a[[paste0(id_var)]],
                                    "_at_",
                                    names(x = at),
                                    "_=_",
                                    a[[names(x = at)]])
      colnames(within_var) <- paste0(a[[paste0(id_var)]],
                                     "_at_",
                                     names(x = at),
                                     "_=_",
                                     a[[names(x = at)]])
      colnames(between_var) <- paste0(a[[paste0(id_var)]],
                                      "_at_",
                                      names(x = at),
                                      "_=_",
                                      a[[names(x = at)]])
      colnames(pooled_se) <- paste0(a[[paste0(id_var)]],
                                    "_at_",
                                    names(x = at),
                                    "_=_",
                                    a[[names(x = at)]])

      for (col in colnames(test$model)) {
        assign(x = paste0("mean_",
                          col),
               value = data.frame(matrix(
                 data = rep(x = NA, times = length(a$fitted)),
                 ncol = length(a$fitted)
               )))

        d <- get(paste0("mean_",
                        col))
        colnames(d) <- paste0(a[[paste0(id_var)]],
                              "_at_",
                              names(x = at),
                              "_=_",
                              a[[names(x = at)]])

        assign(x = paste0("mean_",
                          col),
               value = d)

      }

      for (id in ids) {
        for (value in regpred[[names(at)]]) {
          mean_pred[[paste0(id,
                            "_at_",
                            names(x = at),
                            "_=_",
                            value)]] <-
            mean(x = pred_coef[[paste0(id,
                                       "_at_",
                                       names(x = at),
                                       "_=_",
                                       value)]])
          within_var[[paste0(id,
                             "_at_",
                             names(x = at),
                             "_=_",
                             value)]] <- mean(x = pred_se[[paste0(id,
                                                                  "_at_",
                                                                  names(x = at),
                                                                  "_=_",
                                                                  value)]] ^ 2)

          between_var[[paste0(id,
                              "_at_",
                              names(x = at),
                              "_=_",
                              value)]] <-
            var(x = pred_coef[[paste0(id,
                                      "_at_",
                                      names(x = at),
                                      "_=_",
                                      value)]])

          pooled_se[[paste0(id,
                            "_at_",
                            names(x = at),
                            "_=_",
                            value)]] <-
            sqrt(within_var[[paste0(id,
                                    "_at_",
                                    names(x = at),
                                    "_=_",
                                    value)]] + between_var[[paste0(id,
                                                                   "_at_",
                                                                   names(x = at),
                                                                   "_=_",
                                                                   value)]] + (between_var[[paste0(id,
                                                                                                   "_at_",
                                                                                                   names(x = at),
                                                                                                   "_=_",
                                                                                                   value)]] / imps))

          for (col in colnames(test$model)) {
            d <- get(paste0("mean_",
                            col))

            d[[paste0(id,
                      "_at_",
                      names(x = at),
                      "_=_",
                      value)]] <-
              mean(x = get(paste0(col))[[paste0(id,
                                                "_at_",
                                                names(x = at),
                                                "_=_",
                                                value)]])

            assign(x = paste0("mean_",
                              col),
                   value = d)


          }

        }
      }

      pred <- data.frame(
        pred = rep(x = NA,
                   times = length(a[[names(x = at)]])),
        pred_se = rep(x = NA,
                      times = length(a[[names(x = at)]]))
      )
      pred[[paste0(id_var)]] <- a[[paste0(id_var)]]

      pred[["name"]] <- paste0(a[[paste0(id_var)]],
                               "_at_",
                               names(x = at),
                               "_=_",
                               a[[names(x = at)]])
      pred[["at_level"]] <- a[[names(x = at)]]

      for (col in colnames(test$model)) {
        pred[[paste0(col)]] <- NA



      }

      for (id in ids) {
        for (value in regpred[[names(at)]]) {
          pred$pred[pred[["name"]] == paste0(id,
                                             "_at_",
                                             names(x = at),
                                             "_=_",
                                             value)] <- mean_pred[[paste0(id,
                                                                          "_at_",
                                                                          names(x = at),
                                                                          "_=_",
                                                                          value)]]
          pred$pred_se[pred[["name"]] == paste0(id,
                                                "_at_",
                                                names(x = at),
                                                "_=_",
                                                value)] <-
            pooled_se[[paste0(id,
                              "_at_",
                              names(x = at),
                              "_=_",
                              value)]]

          for (col in colnames(test$model)) {
            pred[[paste0(col)]][pred[["name"]] == paste0(id,
                                                         "_at_",
                                                         names(x = at),
                                                         "_=_",
                                                         value)] <-
              get(paste0("mean_",
                         col))[[paste0(id,
                                       "_at_",
                                       names(x = at),
                                       "_=_",
                                       value)]]

          }
        }
      }

      return(pred)






    }


  }

  ## for hierarchical linear regressions ----

  if (type == "lmer") {

    ifelse(!require(lme4),
           install.packages("lme4", dep = T),
           require(lme4))
    require(lme4)

    ifelse(!require(merTools),
           install.packages("merTools", dep = T),
           require(merTools))
    require(merTools)



    if (is.null(at)) {
      ifelse(!require(prediction),
             install.packages("prediction", dep = T),
             require(prediction))
      require(prediction)

      test  <- lmer(formula = formula,
                    data = df[df[[imp_var]] == 1, ])

      b <- summary(test)

      a <- prediction::prediction(model = test,
                                  at = at,
                                  calculate_se = T,
                                  type = "response")

      N <- length(a$fitted)
      imps <- length(unique(df[[imp_var]]))
      ids <- unique(df[[paste0(id_var)]])

      ses <- predictInterval(merMod = test,
                             newdata = df,
                             seed = bootSeed,
                             n.sims = bootSims,
                             which = "full",
                             stat = "mean",
                             include.resid.var = T)

      ses$se <- (ses$upr - ses$lwr) / (2*1.96)



      for (col in variables) {
        assign(x = paste0(col),
               value = data.frame(matrix(
                 data = rep(x = NA, times = imps * length(a$fitted)),
                 ncol = length(a$fitted)
               )))

        d <- get(paste0(col))
        colnames(d) <- a[[paste0(id_var)]]
        d$imputation <- c(1:imps)

        assign(x = paste0(col),
               value = d)

      }

      pred_coef <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(a$fitted)),
          ncol = length(a$fitted)
        ))
      pred_se <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(a$se.fitted)),
          ncol = length(a$se.fitted)
        ))

      colnames(pred_coef) <- ids
      colnames(pred_se) <- ids

      pred_coef$imputation <- c(1:imps)
      pred_se$imputation <- c(1:imps)

      for (imp in unique(df[[imp_var]])) {
        reg <- lmer(formula = formula,
                    data = df[df[[imp_var]] == imp, ])

        regpred <- prediction::prediction(model = reg,
                                          at = at,
                                          calculate_se = T,
                                          type = "response")

        ses <- predictInterval(merMod = reg,
                               newdata = df[df[[imp_var]] == imp, ],
                               seed = bootSeed,
                               n.sims = bootSims,
                               which = "full",
                               stat = "mean",
                               include.resid.var = T,
                               level = 0.95)

        ses$se <- (ses$upr - ses$lwr) / (2*1.96)

        regmodel <- data.frame(matrix(data = NA,
                                      nrow = length(a$fitted),
                                      ncol = length(variables)))
        colnames(regmodel) <- variables
        regmodel$idvariable <- df[df[[imp_var]] == imp, ][[paste0(id_var)]]
        regmodel$pred <- ses$fit
        regmodel$pred_se <- ses$se

        for (var in variables) {

          regmodel[[paste0(var)]] <- regpred[[paste0(var)]]

        }

        for (id in ids) {
          pred_coef[[paste0(id)]][pred_coef$imputation == imp] <-
            regmodel$pred[regmodel$idvariable == id]
          pred_se[[paste0(id)]][pred_se$imputation == imp] <-
            regmodel$pred_se[regmodel$idvariable == id]

          for (col in variables) {
            d <- get(paste0(col))

            d[[paste0(id)]][get(paste0(col))[["imputation"]] == imp] <-
              regmodel[[paste0(col)]][regmodel$idvariable == id]

            assign(x = paste0(col),
                   value = d)

          }

        }



      }

      # calculate between variance

      mean_pred <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      within_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      between_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      pooled_se <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))

      colnames(mean_pred) <- ids
      colnames(within_var) <- ids
      colnames(between_var) <- ids
      colnames(pooled_se) <- ids

      for (col in variables) {
        assign(x = paste0("mean_",
                          col),
               value = data.frame(matrix(
                 data = rep(x = NA, times = length(a$fitted)),
                 ncol = length(a$fitted)
               )))

        d <- get(paste0("mean_",
                        col))
        colnames(d) <- a[[paste0(id_var)]]

        assign(x = paste0("mean_",
                          col),
               value = d)

      }

      for (id in ids) {
        mean_pred[[paste0(id)]] <- mean(x = pred_coef[[paste0(id)]])
        within_var[[paste0(id)]] <-
          mean(x = pred_se[[paste0(id)]] ^ 2)
        between_var[[paste0(id)]] <-
          var(x = pred_coef[[paste0(id)]])
        pooled_se[[paste0(id)]] <-
          sqrt(within_var[[paste0(id)]] + between_var[[paste0(id)]] + (between_var[[paste0(id)]] / imps))

        for (col in variables) {
          d <- get(paste0("mean_",
                          col))

          d[[paste0(id)]] <-
            mean(x = get(paste0(col))[[paste0(id)]])

          assign(x = paste0("mean_",
                            col),
                 value = d)


        }

      }

      pred <- data.frame(
        pred = rep(x = NA,
                   times = length(ids)),
        pred_se = rep(x = NA,
                      times = length(ids))
      )
      pred[[paste0(id_var)]] <- ids

      for (col in variables) {
        pred[[paste0(col)]] <- NA



      }

      for (id in ids) {
        pred$pred[pred[[paste0(id_var)]] == id] <- mean_pred[[paste0(id)]]
        pred$pred_se[pred[[paste0(id_var)]] == id] <-
          pooled_se[[paste0(id)]]

        for (col in variables) {
          pred[[paste0(col)]][pred[[paste0(id_var)]] == id] <-
            get(paste0("mean_",
                       col))[[paste0(id)]]

        }

      }

      return(pred)




    } else {
      # for specified at() levels


      ifelse(!require(prediction),
             install.packages("prediction", dep = T),
             require(prediction))
      require(prediction)

      test  <- lmer(formula = formula,
                    data = df[df[[imp_var]] == 1, ])

      b <- summary(test)

      a <- prediction::prediction(model = test,
                                  at = at,
                                  calculate_se = T,
                                  type = "response")

      ses <- predictInterval(merMod = test,
                             newdata = df,
                             seed = bootSeed,
                             n.sims = bootSims,
                             which = "full",
                             stat = "mean",
                             include.resid.var = T)


      N <- length(a$fitted)
      imps <- length(unique(df[[imp_var]]))
      ids <- unique(df[[paste0(id_var)]])

      for (col in variables) {
        assign(x = paste0(col),
               value = data.frame(matrix(
                 data = rep(x = NA, times = imps * length(a$fitted)),
                 ncol = length(a$fitted)
               )))

        d <- get(paste0(col))
        colnames(d) <- paste0(ids,
                              "_at_",
                              names(x = at),
                              "_=_",
                              a[[names(x = at)]])
        d$imputation <- c(1:imps)

        assign(x = paste0(col),
               value = d)

      }

      pred_coef <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(a$fitted)),
          ncol = length(a$fitted)
        ))
      pred_se <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(a$se.fitted)),
          ncol = length(a$se.fitted)
        ))

      colnames(pred_coef) <- paste0(ids,
                                    "_at_",
                                    names(x = at),
                                    "_=_",
                                    a[[names(x = at)]])
      colnames(pred_se) <- paste0(ids,
                                  "_at_",
                                  names(x = at),
                                  "_=_",
                                  a[[names(x = at)]])

      pred_coef$imputation <- c(1:imps)
      pred_se$imputation <- c(1:imps)

      for (imp in unique(df[[imp_var]])) {
        reg <- lmer(formula = formula,
                    data = df[df[[imp_var]] == imp, ])

        length_at <- length(at)
        times_at <- 0

        df_imp <- df[df[[imp_var]] == imp, ]

        for (i in 1:length(at)) {

          assign(x = paste0("length_",
                            i),
                 value = length(at[[i]]))

          times_at <- times_at + get(paste0("length_",
                                            i))

          for (j in 1:get(paste0("length_",
                                 i))) {

            d <- df_imp
            d[[names(at)[i]]] <- at[[i]][j]

            assign(x = paste0("df_",
                              i,
                              "_",
                              j),
                   value = d)


          }

        }

        df_se <- df_imp[0,]

        for (i in 1:length(at)) {

          for (j in 1:get(paste0("length_",
                                 i))) {


            df_se <- rbind(df_se,
                           get(paste0("df_",
                                      i,
                                      "_",
                                      j)))


          }

        }




        regpred <- prediction::prediction(model = reg,
                                          at = at,
                                          calculate_se = T,
                                          type = "response")

        ses <- predictInterval(merMod = reg,
                               newdata = df_se,
                               seed = bootSeed,
                               n.sims = bootSims,
                               which = "full",
                               stat = "mean",
                               include.resid.var = T,
                               level = 0.95)

        ses$se <- (ses$upr - ses$lwr) / (2*1.96)

        regmodel <- data.frame(matrix(data = NA,
                                      nrow = length(a$fitted),
                                      ncol = length(variables)))
        colnames(regmodel) <- variables
        regmodel$idvariable <- df[df[[imp_var]] == imp, ][[paste0(id_var)]]
        regmodel$pred <- ses$fit
        regmodel$pred_se <- ses$se
        regmodel$at_level <- regpred[[names(at)]]

        for (var in variables) {

          regmodel[[paste0(var)]] <- regpred[[paste0(var)]]

        }



        for (id in ids) {
          for (value in regpred[[names(at)]]) {
            pred_coef[[paste0(id,
                              "_at_",
                              names(x = at),
                              "_=_",
                              value)]][pred_coef$imputation == imp] <-
              regmodel$pred[(regmodel$idvariable == id &
                               regmodel$at_level == value)]

            pred_se[[paste0(id,
                            "_at_",
                            names(x = at),
                            "_=_",
                            value)]][pred_se$imputation == imp] <-
              regmodel$pred_se[(regmodel$idvariable == id &
                                  regmodel$at_level == value)]

            for (col in variables) {
              d <- get(paste0(col))

              d[[paste0(id,
                        "_at_",
                        names(x = at),
                        "_=_",
                        value)]][get(paste0(col))[["imputation"]] == imp] <-
                regmodel[[paste0(col)]][(regmodel$idvariable == id & regmodel$at_level == value)]

              assign(x = paste0(col),
                     value = d)

            }


          }


        }



      }

      # calculate between variance

      mean_pred <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      within_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      between_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      pooled_se <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))

      colnames(mean_pred) <- paste0(ids,
                                    "_at_",
                                    names(x = at),
                                    "_=_",
                                    a[[names(x = at)]])
      colnames(within_var) <- paste0(ids,
                                     "_at_",
                                     names(x = at),
                                     "_=_",
                                     a[[names(x = at)]])
      colnames(between_var) <- paste0(ids,
                                      "_at_",
                                      names(x = at),
                                      "_=_",
                                      a[[names(x = at)]])
      colnames(pooled_se) <- paste0(ids,
                                    "_at_",
                                    names(x = at),
                                    "_=_",
                                    a[[names(x = at)]])

      for (col in variables) {
        assign(x = paste0("mean_",
                          col),
               value = data.frame(matrix(
                 data = rep(x = NA, times = length(a$fitted)),
                 ncol = length(a$fitted)
               )))

        d <- get(paste0("mean_",
                        col))
        colnames(d) <- paste0(a[[paste0(id_var)]],
                              "_at_",
                              names(x = at),
                              "_=_",
                              a[[names(x = at)]])

        assign(x = paste0("mean_",
                          col),
               value = d)

      }

      for (id in ids) {
        for (value in regpred[[names(at)]]) {
          mean_pred[[paste0(id,
                            "_at_",
                            names(x = at),
                            "_=_",
                            value)]] <-
            mean(x = pred_coef[[paste0(id,
                                       "_at_",
                                       names(x = at),
                                       "_=_",
                                       value)]])
          within_var[[paste0(id,
                             "_at_",
                             names(x = at),
                             "_=_",
                             value)]] <- mean(x = pred_se[[paste0(id,
                                                                  "_at_",
                                                                  names(x = at),
                                                                  "_=_",
                                                                  value)]] ^ 2)

          between_var[[paste0(id,
                              "_at_",
                              names(x = at),
                              "_=_",
                              value)]] <-
            var(x = pred_coef[[paste0(id,
                                      "_at_",
                                      names(x = at),
                                      "_=_",
                                      value)]])

          pooled_se[[paste0(id,
                            "_at_",
                            names(x = at),
                            "_=_",
                            value)]] <-
            sqrt(within_var[[paste0(id,
                                    "_at_",
                                    names(x = at),
                                    "_=_",
                                    value)]] + between_var[[paste0(id,
                                                                   "_at_",
                                                                   names(x = at),
                                                                   "_=_",
                                                                   value)]] + (between_var[[paste0(id,
                                                                                                   "_at_",
                                                                                                   names(x = at),
                                                                                                   "_=_",
                                                                                                   value)]] / imps))

          for (col in variables) {
            d <- get(paste0("mean_",
                            col))

            d[[paste0(id,
                      "_at_",
                      names(x = at),
                      "_=_",
                      value)]] <-
              mean(x = get(paste0(col))[[paste0(id,
                                                "_at_",
                                                names(x = at),
                                                "_=_",
                                                value)]])

            assign(x = paste0("mean_",
                              col),
                   value = d)


          }

        }
      }

      pred <- data.frame(
        pred = rep(x = NA,
                   times = length(a[[names(x = at)]])),
        pred_se = rep(x = NA,
                      times = length(a[[names(x = at)]]))
      )
      pred[[paste0(id_var)]] <- ids

      pred[["name"]] <- paste0(ids,
                               "_at_",
                               names(x = at),
                               "_=_",
                               a[[names(x = at)]])
      pred[["at_level"]] <- a[[names(x = at)]]

      for (col in variables) {
        pred[[paste0(col)]] <- NA



      }

      for (id in ids) {
        for (value in regpred[[names(at)]]) {
          pred$pred[pred[["name"]] == paste0(id,
                                             "_at_",
                                             names(x = at),
                                             "_=_",
                                             value)] <- mean_pred[[paste0(id,
                                                                          "_at_",
                                                                          names(x = at),
                                                                          "_=_",
                                                                          value)]]
          pred$pred_se[pred[["name"]] == paste0(id,
                                                "_at_",
                                                names(x = at),
                                                "_=_",
                                                value)] <-
            pooled_se[[paste0(id,
                              "_at_",
                              names(x = at),
                              "_=_",
                              value)]]

          for (col in variables) {
            pred[[paste0(col)]][pred[["name"]] == paste0(id,
                                                         "_at_",
                                                         names(x = at),
                                                         "_=_",
                                                         value)] <-
              get(paste0("mean_",
                         col))[[paste0(id,
                                       "_at_",
                                       names(x = at),
                                       "_=_",
                                       value)]]

          }
        }
      }

      return(pred)






    }




  }

  # for glm models ----

  if (type == "glm") {
    stopifnot(!is.null(linkfunction))



    if (is.null(at)) {
      ifelse(!require("prediction"),
             install.packages("prediction", dep = T),
             require("prediction"))
      require("prediction")

      test  <- glm(formula = formula,
                   data = df[df[[imp_var]] == 1, ],
                   family = linkfunction)

      a <- prediction::prediction(model = test,
                                  at = at,
                                  calculate_se = T)

      N <- length(a$fitted)
      imps <- length(unique(df[[imp_var]]))
      ids <- unique(df[[paste0(id_var)]])

      for (col in colnames(test$model)) {
        assign(x = paste0(col),
               value = data.frame(matrix(
                 data = rep(x = NA, times = imps * length(a$fitted)),
                 ncol = length(a$fitted)
               )))

        d <- get(paste0(col))
        colnames(d) <- a[[paste0(id_var)]]
        d$imputation <- c(1:imps)

        assign(x = paste0(col),
               value = d)

      }

      pred_coef <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(a$fitted)),
          ncol = length(a$fitted)
        ))
      pred_se <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(a$se.fitted)),
          ncol = length(a$se.fitted)
        ))

      colnames(pred_coef) <- a[[paste0(id_var)]]
      colnames(pred_se) <- a[[paste0(id_var)]]

      pred_coef$imputation <- c(1:imps)
      pred_se$imputation <- c(1:imps)

      for (imp in unique(df[[imp_var]])) {
        reg <- glm(formula = formula,
                   data = df[df[[imp_var]] == imp, ],
                   family = linkfunction)

        regpred <- prediction::prediction(model = reg,
                                          at = at,
                                          calculate_se = T)

        regmodel <- reg$model
        regmodel$idvariable <- regpred[[paste0(id_var)]]

        for (id in ids) {
          pred_coef[[paste0(id)]][pred_coef$imputation == imp] <-
            regpred$fitted[regpred[[paste0(id_var)]] == id]
          pred_se[[paste0(id)]][pred_se$imputation == imp] <-
            regpred$se.fitted[regpred[[paste0(id_var)]] == id]

          for (col in colnames(reg$model)) {
            d <- get(paste0(col))

            d[[paste0(id)]][get(paste0(col))[["imputation"]] == imp] <-
              regmodel[[paste0(col)]][regmodel$idvariable == id]

            assign(x = paste0(col),
                   value = d)

          }

        }



      }

      # calculate between variance

      mean_pred <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      within_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      between_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      pooled_se <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))

      colnames(mean_pred) <- a[[paste0(id_var)]]
      colnames(within_var) <- a[[paste0(id_var)]]
      colnames(between_var) <- a[[paste0(id_var)]]
      colnames(pooled_se) <- a[[paste0(id_var)]]

      for (col in colnames(test$model)) {
        assign(x = paste0("mean_",
                          col),
               value = data.frame(matrix(
                 data = rep(x = NA, times = length(a$fitted)),
                 ncol = length(a$fitted)
               )))

        d <- get(paste0("mean_",
                        col))
        colnames(d) <- a[[paste0(id_var)]]

        assign(x = paste0("mean_",
                          col),
               value = d)

      }

      for (id in ids) {
        mean_pred[[paste0(id)]] <- mean(x = pred_coef[[paste0(id)]])
        within_var[[paste0(id)]] <-
          mean(x = pred_se[[paste0(id)]] ^ 2)
        between_var[[paste0(id)]] <-
          var(x = pred_coef[[paste0(id)]])
        pooled_se[[paste0(id)]] <-
          sqrt(within_var[[paste0(id)]] + between_var[[paste0(id)]] + (between_var[[paste0(id)]] / imps))

        for (col in colnames(test$model)) {
          d <- get(paste0("mean_",
                          col))

          d[[paste0(id)]] <-
            mean(x = get(paste0(col))[[paste0(id)]])

          assign(x = paste0("mean_",
                            col),
                 value = d)


        }

      }

      pred <- data.frame(
        pred = rep(x = NA,
                   times = length(ids)),
        pred_se = rep(x = NA,
                      times = length(ids))
      )
      pred[[paste0(id_var)]] <- ids

      for (col in colnames(test$model)) {
        pred[[paste0(col)]] <- NA



      }

      for (id in ids) {
        pred$pred[pred[[paste0(id_var)]] == id] <- mean_pred[[paste0(id)]]
        pred$pred_se[pred[[paste0(id_var)]] == id] <-
          pooled_se[[paste0(id)]]

        for (col in colnames(test$model)) {
          pred[[paste0(col)]][pred[[paste0(id_var)]] == id] <-
            get(paste0("mean_",
                       col))[[paste0(id)]]

        }

      }

      return(pred)




    } else {
      # for specified at() levels


      ifelse(!require(prediction),
             install.packages("prediction", dep = T),
             require(prediction))
      require(prediction)

      test  <- glm(formula = formula,
                   data = df[df[[imp_var]] == 1, ],
                   family = linkfunction)

      a <- prediction::prediction(model = test,
                                  at = at,
                                  calculate_se = T)

      N <- length(a$fitted)
      imps <- length(unique(df[[imp_var]]))
      ids <- unique(df[[paste0(id_var)]])

      for (col in colnames(test$model)) {
        assign(x = paste0(col),
               value = data.frame(matrix(
                 data = rep(x = NA, times = imps * length(a$fitted)),
                 ncol = length(a$fitted)
               )))

        d <- get(paste0(col))
        colnames(d) <- paste0(a[[paste0(id_var)]],
                              "_at_",
                              names(x = at),
                              "_=_",
                              a[[names(x = at)]])
        d$imputation <- c(1:imps)

        assign(x = paste0(col),
               value = d)

      }

      pred_coef <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(a$fitted)),
          ncol = length(a$fitted)
        ))
      pred_se <-
        data.frame(matrix(
          data = rep(x = NA, times = imps * length(a$se.fitted)),
          ncol = length(a$se.fitted)
        ))

      colnames(pred_coef) <- paste0(a[[paste0(id_var)]],
                                    "_at_",
                                    names(x = at),
                                    "_=_",
                                    a[[names(x = at)]])
      colnames(pred_se) <- paste0(a[[paste0(id_var)]],
                                  "_at_",
                                  names(x = at),
                                  "_=_",
                                  a[[names(x = at)]])

      pred_coef$imputation <- c(1:imps)
      pred_se$imputation <- c(1:imps)

      for (imp in unique(df[[imp_var]])) {
        reg <- glm(formula = formula,
                   data = df[df[[imp_var]] == imp, ],
                   family = linkfunction)

        regpred <- prediction::prediction(model = reg,
                                          at = at,
                                          calculate_se = T)

        regmodel <- reg$model
        regmodel$idvariable <- unique(regpred[[paste0(id_var)]])

        for (id in ids) {
          for (value in regpred[[names(at)]]) {
            pred_coef[[paste0(id,
                              "_at_",
                              names(x = at),
                              "_=_",
                              value)]][pred_coef$imputation == imp] <-
              regpred$fitted[(regpred[[paste0(id_var)]] == id &
                                regpred[[names(at)]] == value)]

            pred_se[[paste0(id,
                            "_at_",
                            names(x = at),
                            "_=_",
                            value)]][pred_se$imputation == imp] <-
              regpred$se.fitted[(regpred[[paste0(id_var)]] == id &
                                   regpred[[names(at)]] == value)]

            for (col in colnames(reg$model)) {
              d <- get(paste0(col))

              d[[paste0(id,
                        "_at_",
                        names(x = at),
                        "_=_",
                        value)]][get(paste0(col))[["imputation"]] == imp] <-
                regmodel[[paste0(col)]][regmodel$idvariable == id]

              assign(x = paste0(col),
                     value = d)

            }


          }


        }



      }

      # calculate between variance

      mean_pred <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      within_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      between_var <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))
      pooled_se <-
        data.frame(matrix(
          data = rep(x = NA, times = length(a$fitted)),
          ncol = length(a$fitted)
        ))

      colnames(mean_pred) <- paste0(a[[paste0(id_var)]],
                                    "_at_",
                                    names(x = at),
                                    "_=_",
                                    a[[names(x = at)]])
      colnames(within_var) <- paste0(a[[paste0(id_var)]],
                                     "_at_",
                                     names(x = at),
                                     "_=_",
                                     a[[names(x = at)]])
      colnames(between_var) <- paste0(a[[paste0(id_var)]],
                                      "_at_",
                                      names(x = at),
                                      "_=_",
                                      a[[names(x = at)]])
      colnames(pooled_se) <- paste0(a[[paste0(id_var)]],
                                    "_at_",
                                    names(x = at),
                                    "_=_",
                                    a[[names(x = at)]])

      for (col in colnames(test$model)) {
        assign(x = paste0("mean_",
                          col),
               value = data.frame(matrix(
                 data = rep(x = NA, times = length(a$fitted)),
                 ncol = length(a$fitted)
               )))

        d <- get(paste0("mean_",
                        col))
        colnames(d) <- paste0(a[[paste0(id_var)]],
                              "_at_",
                              names(x = at),
                              "_=_",
                              a[[names(x = at)]])

        assign(x = paste0("mean_",
                          col),
               value = d)

      }

      for (id in ids) {
        for (value in regpred[[names(at)]]) {
          mean_pred[[paste0(id,
                            "_at_",
                            names(x = at),
                            "_=_",
                            value)]] <-
            mean(x = pred_coef[[paste0(id,
                                       "_at_",
                                       names(x = at),
                                       "_=_",
                                       value)]])
          within_var[[paste0(id,
                             "_at_",
                             names(x = at),
                             "_=_",
                             value)]] <- mean(x = pred_se[[paste0(id,
                                                                  "_at_",
                                                                  names(x = at),
                                                                  "_=_",
                                                                  value)]] ^ 2)

          between_var[[paste0(id,
                              "_at_",
                              names(x = at),
                              "_=_",
                              value)]] <-
            var(x = pred_coef[[paste0(id,
                                      "_at_",
                                      names(x = at),
                                      "_=_",
                                      value)]])

          pooled_se[[paste0(id,
                            "_at_",
                            names(x = at),
                            "_=_",
                            value)]] <-
            sqrt(within_var[[paste0(id,
                                    "_at_",
                                    names(x = at),
                                    "_=_",
                                    value)]] + between_var[[paste0(id,
                                                                   "_at_",
                                                                   names(x = at),
                                                                   "_=_",
                                                                   value)]] + (between_var[[paste0(id,
                                                                                                   "_at_",
                                                                                                   names(x = at),
                                                                                                   "_=_",
                                                                                                   value)]] / imps))

          for (col in colnames(test$model)) {
            d <- get(paste0("mean_",
                            col))

            d[[paste0(id,
                      "_at_",
                      names(x = at),
                      "_=_",
                      value)]] <-
              mean(x = get(paste0(col))[[paste0(id,
                                                "_at_",
                                                names(x = at),
                                                "_=_",
                                                value)]])

            assign(x = paste0("mean_",
                              col),
                   value = d)


          }

        }
      }

      pred <- data.frame(
        pred = rep(x = NA,
                   times = length(a[[names(x = at)]])),
        pred_se = rep(x = NA,
                      times = length(a[[names(x = at)]]))
      )
      pred[[paste0(id_var)]] <- a[[paste0(id_var)]]

      pred[["name"]] <- paste0(a[[paste0(id_var)]],
                               "_at_",
                               names(x = at),
                               "_=_",
                               a[[names(x = at)]])
      pred[["at_level"]] <- a[[names(x = at)]]

      for (col in colnames(test$model)) {
        pred[[paste0(col)]] <- NA



      }

      for (id in ids) {
        for (value in regpred[[names(at)]]) {
          pred$pred[pred[["name"]] == paste0(id,
                                             "_at_",
                                             names(x = at),
                                             "_=_",
                                             value)] <- mean_pred[[paste0(id,
                                                                          "_at_",
                                                                          names(x = at),
                                                                          "_=_",
                                                                          value)]]
          pred$pred_se[pred[["name"]] == paste0(id,
                                                "_at_",
                                                names(x = at),
                                                "_=_",
                                                value)] <-
            pooled_se[[paste0(id,
                              "_at_",
                              names(x = at),
                              "_=_",
                              value)]]

          for (col in colnames(test$model)) {
            pred[[paste0(col)]][pred[["name"]] == paste0(id,
                                                         "_at_",
                                                         names(x = at),
                                                         "_=_",
                                                         value)] <-
              get(paste0("mean_",
                         col))[[paste0(id,
                                       "_at_",
                                       names(x = at),
                                       "_=_",
                                       value)]]

          }
        }
      }

      return(pred)






    }




  }




}
