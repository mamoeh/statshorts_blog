# ----- Consistent random sample queries using cell keys -----
# M. Möhler, Aug. 2025

library(simPop)

# ----- custom functions -----

## helper function: read out a sample based on a cell key

ck_get_s <- function(rkeys, n, ck) {
  
  N  <- length(rkeys)
  s1 <- which(rkeys >= ck)
  n1 <- length(s1)
  
  # the pseudo-sample based on the cell key includes all records with record-keys
  # larger than the cell key; if they are not enough, we append the rest, starting
  # from 0
  if(n1 >= n) {
    s <- s1[1:n]
  } else {
    s <- c((1:N)[1:(n - n1)], s1)
  }
  
  s
}


## function for random sample queries

rs_query <- function(formula, data, sfrac = 1, rkey = NA, ...) {
  
  N <- nrow(data)
  n <- ceiling(sfrac * N)
  
  if(is.na(rkey)) {
    # CASE 1: default RSQ (no cell keys)
    
    # draw sample
    s <- sample(1:N, n)
    s_data <- data[s, ]
    
    # query based on random sample
    q <- xtabs(formula, data = s_data, ...)
    
  } else {
    # CASE 2: RSQ with cell keys
    
    # order by record keys
    data_rk <- data[order(data[, rkey], decreasing = FALSE), ]
    
    # calculate cell keys
    formula_ck <- update.formula(formula, rkey ~ .)
    q_ck <- xtabs(formula_ck, data = data, ...) %% 1
    
    # for each cell in the target table ...
    q <- xtabs(formula, data, ...)
    for(i in seq(q)) {
      
      # ... draw sample according to cell key ...
      s_i <- ck_get_s(data_rk[, rkey], n, q_ck[i])
      
      # ... tabulate statistic from query
      q_i <- xtabs(formula, data[s_i, ], ...)
      q[i] <- q_i[i]
    }
  }
  
  # multiply by sample weight
  round(q * (1/sfrac))
}


# ----- application -----

## (1) prepare data

data("eusilcS")

# keep only records with age and economic status information
eusilcS <- eusilcS[!is.na(eusilcS$age) & !is.na(eusilcS$pl030), ]
# group age
eusilcS$age <- cut(eusilcS$age, c(0, 18, 30, 65, 80, 100))
# group economic status
econ_labels <- c("full time or part time", "unemployed or in training", "retired, inactive, or in care work")
eusilcS$econ_stat <- cut(as.numeric(eusilcS$pl030), c(0, 2, 4, 7), right = TRUE, labels = econ_labels)
# append record keys
set.seed(42)
eusilcS$rkey <- runif(nrow(eusilcS))


## (2) run queries

# full query
xtabs(~ econ_stat + age, data = eusilcS)
rs_query(~ econ_stat + age, data = eusilcS, sfrac = 1)
# marginal query
xtabs(~ age,       data = eusilcS)
xtabs(~ econ_stat, data = eusilcS)

# RS query showing inconsistency over multiple iterations
set.seed(20250824)
rs_query(~ econ_stat + age, data = eusilcS, sfrac = 0.8)
rs_query(~ econ_stat + age, data = eusilcS, sfrac = 0.8)
# RS query showing inconsistency over inner and marginal queries
rs_query(~ age,       data = eusilcS, sfrac = 0.8)
rs_query(~ econ_stat, data = eusilcS, sfrac = 0.8)

# cell keys
(ck_inner <- round(xtabs(rkey ~ econ_stat + age, data = eusilcS), 6) %% 1)
(ck_margn <- round(xtabs(rkey ~ age,             data = eusilcS), 6) %% 1)
# starting IDs of RSQs
rkeys <- sort(eusilcS$rkey, decreasing = FALSE)
mapply(ck_inner, FUN = function(x){min(which(rkeys > x))})
mapply(ck_margn, FUN = function(x){min(which(rkeys > x))})

# RS query made consistent by using cell keys
# count queries
(rs_ct_inner <- rs_query(~ econ_stat + age, data = eusilcS, sfrac = 0.8, rkey = "rkey"))
(rs_ct_margn <- rs_query(~ age,             data = eusilcS, sfrac = 0.8, rkey = "rkey"))
# sum queries (net income)
(rs_sm_inner <- rs_query(netIncome ~ econ_stat + age, data = eusilcS, sfrac = 0.8, rkey = "rkey"))
(rs_sm_margn <- rs_query(netIncome ~ age,             data = eusilcS, sfrac = 0.8, rkey = "rkey"))
# calculated mean (net income)
round(rs_sm_inner / rs_ct_inner)
round(rs_sm_margn / rs_ct_margn)

