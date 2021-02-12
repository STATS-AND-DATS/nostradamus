#' A Predictor Function
#'
#' This function allows you to find the best predictor
#' @keywords predictor
#' @export
#' @import rvest, dplyr, stringr, tidyr
#' @examples
#' best_predictor(dataframe, response)


best_predictor <- function(dataframe, response) {
  if (sum(sapply(dataframe, function(x) {is.numeric(x) | is.factor(x)})) < ncol(dataframe)) {
    stop("Make sure that all variables are of class numeric/factor!")
  }
  # pre-allocate vectors
  varname <- c()
  vartype <- c()
  R2 <- c()
  R2_log <- c()
  R2_quad <- c()
  AIC <- c()
  AIC_log <- c()
  AIC_quad <- c()
  y <- dataframe[ ,response]

  # # # # # NUMERIC RESPONSE # # # # #
  if (is.numeric(y)) {

    for (i in 1:ncol(dataframe)) {

      x <- dataframe[ ,i]
      varname[i] <- names(dataframe)[i]

      if (class(x) %in% c("numeric", "integer")) {
        vartype[i] <- "numeric"
      } else {
        vartype[i] <- "categorical"
      }

      if (!identical(y, x)) {

        # linear: y ~ x
        R2[i] <- summary(lm(y ~ x))$r.squared

        # log-transform: y ~ log(x)
        if (is.numeric(x)) {
          if (min(x) <= 0) { # if y ~ log(x) for min(x) <= 0, do y ~ log(x + abs(min(x)) + 1)
            R2_log[i] <- summary(lm(y ~ log(x + abs(min(x)) + 1)))$r.squared
          } else {
            R2_log[i] <- summary(lm(y ~ log(x)))$r.squared
          }
        } else {
          R2_log[i] <- NA
        }

        # quadratic: y ~ x + x^2
        if (is.numeric(x)) {
          R2_quad[i] <- summary(lm(y ~ x + I(x^2)))$r.squared
        } else {
          R2_quad[i] <- NA
        }

      } else {
        R2[i] <- NA
        R2_log[i] <- NA
        R2_quad[i] <- NA
      }
    }

    print(paste("Response variable:", response))

    data.frame(varname,
               vartype,
               R2 = round(R2, 3),
               R2_log = round(R2_log, 3),
               R2_quad = round(R2_quad, 3)) %>%
      mutate(max_R2 = pmax(R2, R2_log, R2_quad, na.rm = T)) %>%
      arrange(desc(max_R2))


    # # # # # CATEGORICAL RESPONSE # # # # #
  } else {

    for (i in 1:ncol(dataframe)) {

      x <- dataframe[ ,i]
      varname[i] <- names(dataframe)[i]

      if (class(x) %in% c("numeric", "integer")) {
        vartype[i] <- "numeric"
      } else {
        vartype[i] <- "categorical"
      }

      if (!identical(y, x)) {

        # linear: y ~ x
        AIC[i] <- summary(glm(y ~ x, family = "binomial"))$aic

        # log-transform: y ~ log(x)
        if (is.numeric(x)) {
          if (min(x) <= 0) { # if y ~ log(x) for min(x) <= 0, do y ~ log(x + abs(min(x)) + 1)
            AIC_log[i] <- summary(glm(y ~ log(x + abs(min(x)) + 1), family = "binomial"))$aic
          } else {
            AIC_log[i] <- summary(glm(y ~ log(x), family = "binomial"))$aic
          }
        } else {
          AIC_log[i] <- NA
        }

        # quadratic: y ~ x + x^2
        if (is.numeric(x)) {
          AIC_quad[i] <- summary(glm(y ~ x + I(x^2), family = "binomial"))$aic
        } else {
          AIC_quad[i] <- NA
        }

      } else {
        AIC[i] <- NA
        AIC_log[i] <- NA
        AIC_quad[i] <- NA
      }
    }

    print(paste("Response variable:", response))

    data.frame(varname,
               vartype,
               AIC = round(AIC, 3),
               AIC_log = round(AIC_log, 3),
               AIC_quad = round(AIC_quad, 3)) %>%
      mutate(min_AIC = pmin(AIC, AIC_log, AIC_quad, na.rm = T)) %>%
      arrange(min_AIC)
  }
}

determine_predictor <- function(dataframe, response) {
  if (is.numeric(dataframe[,response])) {
    quant_predictors <- best_predictor(dataframe, response) %>%
      select(-c(vartype, max_R2)) %>%
      gather(key = "key", value = "R2", -varname) %>%
      filter(!is.na(R2))

    quant_predictors_order <- best_predictor(dataframe, response) %>%
      select(varname, max_R2) %>%
      filter(!is.na(max_R2)) %>%
      arrange(max_R2) %>%
      pull(varname)

    quant_predictors$varname <- factor(quant_predictors$varname,
                                       ordered = T,
                                       levels = quant_predictors_order)
    result <- quant_predictors
    return(result)
  } else {
    qual_predictors <- best_predictor(dataframe, response) %>%
      select(-c(vartype, min_AIC)) %>%
      gather(key = "key", value = "AIC", -varname) %>%
      filter(!is.na(AIC))

    qual_predictors_order <- best_predictor(dataframe, response) %>%
      select(varname, min_AIC) %>%
      filter(!is.na(min_AIC)) %>%
      arrange(desc(min_AIC)) %>%
      pull(varname)

    qual_predictors$varname <- factor(qual_predictors$varname,
                                      ordered = T,
                                      levels = qual_predictors_order)
  }
  result <- qual_predictors
  return(result)

}
