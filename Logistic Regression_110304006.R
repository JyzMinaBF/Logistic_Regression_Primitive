# This is the first function.
# It is used to automatically generate simulated data and directly do the logistic regression
# to find the beta
# There are six arguments that we can control 
# 1. dim: the number of dimensions that we want to generate for X
# 2. num_samples: number of samples that we want to generate for X
# 3. restr: if the distance between calculated beta n and calculated beta n-1 is lower than this, stop the iteration
# 4. max_times: the maximum times to iterate
# 5. true_beta: the true beta that we use to generate y, and it will be reduced due to dim
# 6. whether to report the final result
logistic_reg_simu = function(dim, num_samples, restr = 0.00000001, 
                             max_times = 100000, true_beta = c(2, 3, 1, 5, 4, 1, 7, 8, 2, 9, 6), 
                             report = F){
  library(MASS)
  
  # the function to implement newton-raphson method 
  newton_raphson = function(x, y, beta){
    y_hat = 1 / (1 + exp(-c(x %*% beta))) #use the old beta to predict first
    first_dev = t(x) %*% (y - y_hat) # get the first order derivation
    second_dev = t(x) %*% diag(y_hat) %*% x # get the second order derivation
    # The code below shows that the new beta is the sum of old beta and first_dev/second_dev
    # (we use solve() to get the inverse to show the division)
    beta_next = beta + solve(second_dev) %*% first_dev 
    return(beta_next)
  }
  
  # the function to calculate the euclidean distance between two beta
  euc_dist = function(x, y){
    if(length(x) == length(y)){ # check
      n = length(x)
      diff = 0
      for(i in 1:n){
        diff = diff + (x[i]-y[i])^2
      }
      dist = diff^(1/2)
      return(dist)
    }
  }
  
  # decide the true beta by dim
  true_beta_here = true_beta[1:(dim+1)]
  # set a difference value that will surely exceed the restriction below
  dist = 1 
  # initialize how many time we iterate
  iterative_times = 0 
  # initialize the value of beta
  test_beta = rep(1, dim+1)  
  #use mvrnorm() to generate X first
  var_mat = matrix(rep(0, dim ^ 2), ncol = dim)
  for(i in 1:dim){
    var_mat[i, i] = 1
  }
  x_inception = rep(1, num_samples) 
  x_samples = cbind(x_inception, mvrnorm(num_samples, rep(0, dim), var_mat))
  # calculate with the true beta to find the probability
  xb = x_samples %*% true_beta_here
  p <- 1/(1 + exp(-xb))
  # use rbinom() with the probability above to generate y 
  y_samples = rbinom(n = num_samples, size = 1, prob = p)
  
  start_time =  as.POSIXct(Sys.time()) # calculate how much time used
  
  # start to iterate to calculate the find the new beta
  while(dist > restr && iterative_times < max_times){
    beta_temp = test_beta # store the value of old beta first
    test_beta = newton_raphson(x_samples, y_samples, test_beta) # get the new beta from the function above
    iterative_times = iterative_times + 1
    dist = euc_dist(test_beta, beta_temp)
  }
  
  # calculate how much time used
  end_time =  as.POSIXct(Sys.time()) 
  used_time = as.numeric(difftime(end_time, start_time, units = "secs")) 
  
  # use the control of argument "report" to decide whether we should generate the result
  if(report == T){
    cat("The value of true beta:", true_beta_here, "\n")
    cat("The value of calculated beta:", test_beta, "\n")
    cat("Itearative Times:", iterative_times, "\n")
    cat("Running Time:", used_time, "s\n")
  }

  # output of the function
  reg_result = list(x_samples = x_samples, y_samples = y_samples, 
                    test_beta = test_beta, 
                    iterative_times = iterative_times, used_time = used_time)
  return(reg_result)
}

# This is the second function.
# It is used to generate the ROC curve of a single simulation of logistic regression
# There are three arguments that we can control
# 1. dim: the same as what we see above 
# 2. num_samples: the same as what we see above
# 3. cal_auc: whether to calculate AUC and show it at the end
roc_curve = function(dim, num_samples, cal_auc = F){
  exper_result = logistic_reg_simu(dim = dim, num_samples = num_samples)
  y_pred =  1/(1 + exp(-(exper_result$x_samples %*% exper_result$test_beta)))
  y_true = exper_result$y_samples
  thre_vector = seq(0, 1, length = 5000) 
  roc_points = data.frame(matrix(c(0, 0), nrow = 1))
  for(thre in thre_vector){
    y_pred_trans = ifelse(y_pred >= thre, 1, 0)
    conf_mat = matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)
    for(i in 1:length(y_pred_trans)){
      if(y_pred_trans[i] == 1 && y_true[i] == 1){
        conf_mat[1, 1] = conf_mat[1, 1] + 1
      }else if(y_pred_trans[i] == 1 && y_true[i] == 0){
        conf_mat[1, 2] = conf_mat[1, 2] + 1
      }else if(y_pred_trans[i] == 0 && y_true[i] == 1){
        conf_mat[2, 1] = conf_mat[2, 1] + 1
      }else if(y_pred_trans[i] == 0 && y_true[i] == 0){
        conf_mat[2, 2] = conf_mat[2, 2] + 1
      }
    }
    tpr = conf_mat[1, 1]/(conf_mat[1 ,1]+conf_mat[2 ,1]) 
    fpr = conf_mat[1, 2]/(conf_mat[1 ,2]+conf_mat[2 ,2])
    roc_points[nrow(roc_points)+1,] = c(fpr, tpr)
  }
  roc_points = roc_points[-1,]
  plot(roc_points[,1], roc_points[,2], 
       type = "l", lwd = 3, col = "red", 
       xlab = "FPR", ylab = "TPR", 
       xlim = c(0.03, 0.97), ylim = c(0.03, 0.97))
  title(main = paste("ROC Curve of Logistic Regression \n with", 
                     as.character(dim),"Predicted Variables and",
                     as.character(num_samples), "Samples"))
  if(cal_auc == T){
    auc = 0
    for(i in 1:4999){
      auc = auc + (roc_points[i,2]+roc_points[i+1,2])*(roc_points[i, 1]-roc_points[i+1, 1])/2
    }
    cat("AUC:", auc, "\n")
  }
}

# This is the third function
# It is used to repeatedly run logistic regression with different dim and num_samples to estimate the difference
# There are four arguments that we can control
# 1. test_times: how many times to repeatedly test for every pair of number of dimensions and number of samples
# 2. test_dim: the numbers of dimensions that you are willing to test
# 3. test_num_samples: the numbers of samples that you are willing to test
# 4. whether to report that it is running on which pair for which times, it can use to realize the speed and estimate when to finish
rep_log_reg_sim = function(test_times = 100, 
                           test_dim = c(2, 3, 4), test_num_samples = c(100, 300), 
                           cond_report = F){
  result = list()
  for (dim in test_dim){
    for (num_samples in test_num_samples){
      
      df_test = data.frame(t(data.frame(rep(0, dim+3))))
      cname = c()
      for(i in 1:(dim+1)){
        cname = c(cname, paste("beta", as.character(i-1), sep = ""))
      }
      colnames(df_test) = c(cname, "Iterative_times", "Running_Time")
      
      for (n in 1:test_times){
        reg_result = logistic_reg_simu(dim, num_samples)

        df_test[nrow(df_test) + 1, ] = c(reg_result$test_beta, reg_result$iterative_times, reg_result$used_time)
        
        if(cond_report == T){
          cat(dim, num_samples, n, "\n")
        }
      }
      df_test = df_test[-1,]
      result[[paste(as.character(dim), ", ", as.character(num_samples), sep = "")]] = df_test
    }
  } 
  return(result)
}

# below is the main working part
true_beta = c(2, 3, 1, 5, 4, 1, 7, 8, 2, 9, 6)
one_time_result = logistic_reg_simu(dim = 2, num_samples = 100, report = T, true_beta = c(2, 3, 1, 5, 4, 1, 7, 8, 2, 9, 6))
roc_curve(dim = 2, num_samples = 2000, cal_auc = T)
result = rep_log_reg_sim(test_dim = c(2, 3, 4), test_num_sampless = c(100, 300), cond_report = T)

# this one is used to show the summary of repeatedly running logistic regression with different dim and num_samples
# we can even see the iterative times and running time is this summary report
for (dim in c(2, 3, 4)){
  for (num_samples in c(100, 300)){
    cat("Dim:", dim, "\n")
    cat("Number of Samples:", num_samples, "\n")
    cat("True beta:", true_beta[1:(dim+1)], "\n")
    res_df = result[[paste(as.character(dim), ", ", as.character(num_samples), sep = "")]]
    for(col in colnames(res_df)){
      cat(col, ":\n","Mean:", mean(c(res_df[,col])), "Var:", var(c(res_df[,col])), "\n")
    }
    cat("\n")
  }
}



