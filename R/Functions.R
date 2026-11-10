###########################################################################################
######################################## FUNCTIONS ########################################
###########################################################################################


# Function to install and/or run the necessary packages -----------------------------------

package.check <- function(packages) {
  packages %in% rownames(installed.packages()) ->
    installed_libraries
  
  if (any(installed_libraries == F)) {
    install.packages(packages[!installed_libraries])
  }
  
  invisible(lapply(packages, library,
                   character.only = T))
}

saveRDS(package.check, "Functions/package.check.rds")


# Function to carried out the analysis of Cruz, Torres & Vencovsky ------------------------

ctv <-
  function(data,
           genotype = "CLONE",
           env = "ENV",
           block = "REP",
           y = "ATR") {
    
    # Transforma em fator
    data1 <- data %>%
      as_tibble() %>%
      mutate_at(c(genotype, env, block), as.factor)
    
    data1 |> 
      group_by(!!sym(genotype)) |> 
      summarise_at(vars(all_of(y)), mean) -> mean_var
    
    mean_var |>
      select_at(vars(all_of(genotype))) |>
      pull() -> names_gen
    
    mean_var |>
      select_at(vars(all_of(y))) |>
      pull() -> mean_gen
    
    names(mean_gen) <- names_gen
    
    #b0 <- mean(data1[[y]])
    
    # ANOVA (desconsiderar os testes)
    fit_anova <- map(paste(y, "~", genotype, "*", env, "+", block, ":", env),
                     ~ anova(lm(., data1)))[[1]] %>% broom::tidy()
    
    # Soma de quadrado e Quadrado médio do resíduo
    SQR <- as_vector(fit_anova[5, 3])
    QMR <- as_vector(fit_anova[5, 4])
    
    
    # Número de repetições
    r <- max(as.numeric(data1[[block]]))
    g1 <- ensym(genotype)
    g2 <- ensym(env)
    g3 <- ensym(y)
    
    # Dados de médias (ainda desordenados)
    data2 <- data %>%
      as_tibble() %>%
      dplyr::select(c(all_of(genotype), all_of(env), all_of(y))) %>%
      group_by(!!g1,!!g2) %>%
      summarise_all(mean) %>%
      pivot_wider(names_from = !!g2, values_from = !!g3) %>%
      ungroup() %>%
      mutate_at(genotype, as.factor)
    
    # Obtendo o Índice Ambiental
    I_unord <- colMeans(data2[, -1]) - mean(as.matrix(data2[, -1]))
    I_ord <- sort(I_unord)
    
    fav <- function(x) {
      if (x < 0) {
        f  <- "Favorable"
      } else {
        f <- "Unfavorable"
      }
      return(f)
    }
    
    Type <- sapply(I_unord, fav)
    
    
    # Removendo a coluna de grupos e organizando as colunas pela ordem dos índices ambientais
    g <- data2[, 1]
    
    env_names <- names(data2)[-1][order(I_unord)]
    
    data3 <- data2[, -1][, order(I_unord)]
    
    
    data4 <- data2[, -1]
    
    
    b0 <- mean(as.matrix(data4))
    mean_geral <- round(rowMeans(data4), 2)
    mean_unfav <- round(rowMeans(data4[, I_unord < 0]), 2)
    mean_fav <- round(rowMeans(data4[, I_unord >= 0]), 2)
    
    
    n_environ <- ncol(data4)
    n_genotype <- nrow(data4)
    n_total <- n_genotype * n_environ * r
    
    gl_b_a <- n_environ * (r - 1)       # gl bloco dentro de ambiente
    gl_a <- n_environ - 1               # gl ambiente
    gl_g <- n_genotype - 1              # gl genótipos
    gl_ga <- gl_a * gl_g                # gl interação genótipo ambiente
    gl_ag <- gl_a * n_genotype          # gl ambiente dentro de genótipo
    gl_ag_i <- gl_a                     # gl ambiente dentro de cada genótipo
    gl <- n_environ * (r - 1) * gl_g    # gl do resíduo
    gl_dev_i <- n_environ - 3
    
    Yi. <- rowSums(data4)              # totais de genótipos
    Y.j <- colSums(data4)              # totais de ambiente
    total <- sum(data4)                # total geral
    
    CC <- ((total * r) ^ 2) / n_total
    
    
    SQAG_i <-
      NULL # Soma de quadrados ambiente dentro de cada genótipo
    for (i in 1:nrow(data4)) {
      x <- as_vector(data4[i, ])
      SQAG_i[i] <- r * (sum(x ^ 2) - (sum(x) ^ 2) / ncol(data4))
    }
    
    
    QMAG_i <-
      SQAG_i / gl_ag_i # Quadrado Ambiente dentro de cada genótipo
    
    F_QMAG_i <- QMAG_i / QMR
    
    SQAG <-
      sum(SQAG_i) # Soma de quadrados Ambiente dentro de genótipo
    
    
    SQGA_j <- NULL # Soma de quadrados genótipo dentro cada ambiente
    for (j in 1:ncol(data4)) {
      x <- as_vector(data4[, j])
      SQGA_j[j] <- r * (sum(x ^ 2) - (sum(x) ^ 2) / nrow(data4))
    }
    
    QMGA_j <-
      SQGA_j / gl_g # Quadrado médio Genótipo dentro de cada ambiente
    
    
    TIj_func <- function(x) {
      if (x < 0) {
        x = 0
      } else {
        x = x - mean(I_unord[I_unord > 0])
      }
      return(x)
    }
    
    TIj <- sapply(I_unord, TIj_func)
    
    
    SQReg_i <- NULL
    for (i in 1:nrow(data4)) {
      x <- as_vector(data4[i, ])
      p1 <-
        (r * (t(x) %*% I_unord - t(x) %*% TIj) ^ 2) / (t(I_unord) %*% I_unord - t(TIj) %*%
                                                         TIj)
      p2 <- (r * (t(x) %*% TIj) ^ 2) / (t(TIj) %*% TIj)
      SQReg_i[i] <- p1 + p2
    }
    
    QMReg_i <- SQReg_i / 2 # Quadrados médios de cada regressão
    
    SQDesv_i <-
      SQAG_i - SQReg_i # Soma de quadrados dos desvios de cada regressão
    
    QMDesv_i <-
      SQDesv_i / (n_environ - 3) # Quadrados Médios dos desvios
    
    desvio_i <- (QMDesv_i - QMR) / r # Estimativas dos desvios
    
    F_dev <- QMDesv_i / QMR # Estatística F para os desvios
    
    prob_dev <- round(pf(F_dev, gl_dev_i, gl, lower.tail = F) * 100, 4)
    
    R2 <- round(SQReg_i / SQAG_i * 100, 2)
    
    X <- cbind(1, I_unord, TIj)
    XtX <- t(X) %*% X
    inv_XtX <- solve(XtX)
    
    # Estimativas dos betas
    b0_b1_b2 <-
      round(inv_XtX %*% t(X) %*% t(as.matrix(data4)), 7) # estima b0, b1 e b2
    
    b1_plus_b2 <-
      as.matrix(data4) %*% TIj / c(t(TIj) %*% TIj) # estima b1+b2
    
    estimates <- rbind(b0_b1_b2, t(b1_plus_b2))
    dimnames(estimates) <-
      list(Betas = c("B0", "B1", "B2", "B1 + B2"),
           Genótipos = row.names(data2))
    
    # Variâncias residual e do desvio da regressão
    s2e <- QMR / r
    
    # Variâncias dos betas
    M_Var <- as.vector(c(diag(inv_XtX * s2e), s2e / t(TIj) %*% TIj))
    names(M_Var) <- c("B0", "B1", "B2", "B1 + B2")
    
    
    t_b1_cal <- as.vector((estimates[2,] - 1) / sqrt(M_Var[[2]]))
    t_b1_b2_cal <- as.vector((estimates[4, ] - 1) / sqrt(M_Var[[4]]))
    
    
    prob_t <- function(x) {
      if (x < 0) {
        prob <- round(200 * (pt(x, gl, lower.tail = T)), 2)
      } else {
        prob <- round(200 * (pt(x, gl, lower.tail = F)), 2)
      }
      return(prob)
    }
    
    prob_b1 <- sapply(t_b1_cal, prob_t)
    
    prob_b2_b2 <- sapply(t_b1_b2_cal, prob_t)
    
    r0 <- tibble(g, data3, S2d = desvio_i) # mudei aqui adicionando g
    
    r1 <-
      data.frame(
        Environment = colnames(data4),
        Mean = colMeans(data2[, -1]),
        Ij = I_unord,
        TIj,
        Type
      )
    
    r2 <- tibble(data2[, 1], mean_geral, mean_unfav, mean_fav)
    colnames(r2) <-
      c(names(r2[, 1]), "Both", "Unfavorable", "Favorable")
    
    r3 <- data.frame(data2[, 1], t(estimates), s2d = desvio_i, R2,  
                     t_b1_cal, t_b1_b2_cal, F_dev,
                     prob_b1, prob_b2_b2, prob_dev) |>
      as_tibble()
    
    
    range_b1 <- range(r3$B1)
    
    range_b1_b2 <- range(r3$B1...B2)
    
    r4 <-
      list(
        b0 = mean_gen,
        df = gl,
        msq = as.vector(QMR),
        r = r,
        I_ord = as.vector(I_ord),
        env_names = env_names,
        genotype = genotype
      )
    
    
    return(
      list(
        data_for_classification = r0,
        mean_gen_per_env = data2,
        mean_env_index = r1,
        mean_gen_in_both_unfav_fav = r2,
        curve_parameters = r3,
        range_b1 = range_b1,
        range_b1_b2 = range_b1_b2,
        beta_variances = round(M_Var, 6),
        parameters_for_simulation = r4
      )
    )
    
  }


saveRDS(ctv, "Functions/ctv.rds")


# Function for data simulation ------------------------------------------------------------


#' @param parameters a list with the following parameters: b0, df, msq, r, I_ord, env_names and, genotypes.
#' @param n_training number of simulated genotypes (default = 500).
#' @param alpha the significance level the estimates of betas and s2d will be tested.

train.set_ctv <- function(parameters = list(),
                          n_training = 500,
                          alpha = 0.05) {
  
    if (length(parameters) == 7 &
        is.numeric(n_training) == T &
        is.numeric(alpha) == T &
        alpha <= 0.1 &
        alpha >= 0.001) {
      
      b0 <- parameters[[1]]
      gl <- parameters[[2]]
      QMR <- parameters[[3]]
      r <- parameters[[4]]
      I_ord <- parameters[[5]]
      env_names <- parameters[[6]]
      genotype <- parameters[[7]]
      
      n_environ <- length(I_ord)
      
      simula_ctv <-
        function(QMR,
                 n_environ,
                 b0,
                 I_b1,
                 I_b1_b2,
                 I_ord) {
          if (I_b1 == 1)
          {
            b1 <- runif(1, 0, 0.89)
          }
          if (I_b1 == 2)
          {
            b1 <- runif(1, 0.90, 1.10)
          }
          if (I_b1 == 3)
          {
            b1 <- runif(1, 1.11, 2)
          }
          
          if (I_b1_b2 == 1)
          {
            b1_b2 <- runif(1, 0, 0.89)
          }
          if (I_b1_b2 == 2)
          {
            b1_b2 <- runif(1, 0.90, 1.10)
          }
          if (I_b1_b2 == 3)
          {
            b1_b2 <- runif(1, 1.11, 2)
          }
          
          b2 = b1_b2 - b1
          
          TIj_func <- function(x) {
            if (x < 0) {
              x = 0
            } else {
              x = x - mean(I_ord[I_ord > 0])
            }
            return(x)
          }
          
          b0_sim <- runif(1, min(b0), max(b0))
          
          TIj <- sapply(I_ord, TIj_func)
          
          erro <- rnorm(n = n_environ,
                        mean = 0,
                        sd = sqrt(QMR))
          
          y_obs <- as.vector(b0_sim + b1 * I_ord + b2 * TIj + erro)
          
          return(list(y_obs = y_obs, TIj = TIj))
        }
      
      # Valores de t e F tabelados acessados dentro de cada repeat
      t_tab.1 = qt(alpha / 2, gl , lower.tail = T)
      t_tab.2 = qt(1 - alpha / 2, gl , lower.tail = T)
      F_tab = qf(1 - alpha, n_environ - 3, gl, lower.tail = T)
      
      # Fun??o para estimar os valores de t
      ts <- function (y_obs, TIj, I_ord, n_environ, QMR, r) {
        # Matrizes
        X <- cbind(1, I_ord, TIj)
        XtX <- t(X) %*% X
        inv_XtX <- solve(XtX)
        
        # Estimativas dos betas
        b0_b1_b2 <-
          as.vector(inv_XtX %*% t(X) %*% y_obs) # estima b0, b1 e b2
        b1_plus_b2 <- t(y_obs) %*% TIj / t(TIj) %*% TIj # estima b1+b2
        
        # Somas de quadrados e quadrados m?dios
        SQAG <- r * (t(y_obs) %*% y_obs - (1 / n_environ) * sum(y_obs) ^ 2)
        SQD <- r * (t(y_obs) %*% y_obs - t(b0_b1_b2) %*% t(X) %*% y_obs)
        SQReg <- SQAG - SQD
        QMD <- SQD / (n_environ - 3)
        
        # R2
        R2 <- SQReg / SQAG
        
        # Vari?ncias residual e do desvio da regress?o
        s2e <- QMR / r
        s2d <- (QMD - QMR) / r
        
        # Vari?ncias dos betas
        M_Var <- inv_XtX * s2e # de b0, b1 e b2
        round(M_Var, 2)
        V_B1_B2 <- QMR / (r * t(TIj) %*% TIj) # de b1 + b2
        
        # t calculado para b1 e b1+b2
        t_b1_cal <- (b0_b1_b2[2] - 1) / sqrt(M_Var[2, 2])
        t_b1_b2_cal <- (b1_plus_b2 - 1) / sqrt(V_B1_B2)
        
        # F calculado
        F_cal <- QMD / QMR
        
        return(
          list(
            t_b1_cal = t_b1_cal,
            t_b1_b2_cal = t_b1_b2_cal,
            F_cal = F_cal,
            s2d = s2d,
            b1 = b0_b1_b2[[2]],
            b1_plus_b2 = b1_plus_b2[[1]]
          )
        )
        
      }
      
      # Criando as condi??es
      u1 <- quote(ts_c[[1]] < t_tab.1)
      u2 <- quote(ts_c[[1]] > t_tab.1)
      u3 <- quote(ts_c[[1]] < t_tab.2)
      u4 <- quote(ts_c[[1]] > t_tab.2)
      
      u5 <- quote(ts_c[[2]] < t_tab.1)
      u6 <- quote(ts_c[[2]] > t_tab.1)
      u7 <- quote(ts_c[[2]] < t_tab.2)
      u8 <- quote(ts_c[[2]] > t_tab.2)
      
      u9 <- substitute(ts_c[[3]] < F_tab)
      u10 <- substitute(ts_c[[3]] > F_tab)
      
      conditions <- list(
        substitute(u1 & u5),
        substitute(u1 & u6 & u7),
        substitute(u1 & u8),
        substitute(u2 & u3 & u5),
        substitute(u2 & u3 & u6 & u7),
        substitute(u2 & u3 & u8),
        substitute(u4 & u5),
        substitute(u4 & u6 & u7),
        substitute(u4 & u8)
      )
      
      
      I_b <- list(c(1, 1),
                  c(1, 2),
                  c(1, 3),
                  c(2, 1),
                  c(2, 2),
                  c(2, 3),
                  c(3, 1),
                  c(3, 2),
                  c(3, 3))
      
      
      include_sim <-
        function(conditions,
                 n_training,
                 n_environ,
                 QMR,
                 b0,
                 I_b,
                 I_ord,
                 g,
                 r,
                 stab) {
          cont <- 0
          d <- matrix(NA, n_training, n_environ + 4, byrow = TRUE)
          i <- 1
          
          if (stab == 1) {
            F_cond = u9 # High stability
          } else {
            F_cond = u10 # Low stability
          }
          
          while (cont < n_training)
          {
            # Simulate the observarions
            z <-
              simula_ctv(QMR,
                         n_environ,
                         b0,
                         I_b1 = I_b[1],
                         I_b1_b2 = I_b[2],
                         I_ord)
            y_obs <- z[[1]]
            TIj <- z[[2]]
            
            # Valores de t calculados
            ts_c <- ts(y_obs, TIj, I_ord, n_environ, QMR, r)
            
            # Conditions to include the observation in the data frame
            if (eval(conditions) & eval(F_cond))
            {
              d[i, ] <- c(g, y_obs, ts_c[[4]], ts_c[[5]], ts_c[[6]])
              cont <- cont + 1
              i <- i + 1
            }
            
          }
          return(d)
        }
      
      data <- NULL
      
      a <- Sys.time()
      for (g in 1:9) {
        d <-
          include_sim(conditions[[g]],
                      n_training,
                      n_environ,
                      QMR,
                      b0,
                      I_b[[g]],
                      I_ord,
                      g,
                      r,
                      stab = 1)
        data <- abind::abind(data, d, along = 1)
      }
      
      for (g in 1:9) {
        d <-
          include_sim(conditions[[g]],
                      n_training,
                      n_environ,
                      QMR,
                      b0,
                      I_b[[g]],
                      I_ord,
                      g + 9,
                      r,
                      stab = 2)
        data <- abind::abind(data, d, along = 1)
      }
      
      b <- Sys.time()
      
      time <- as.numeric(b - a, units = "secs")
      
      
      # Conjunto de dados para treinamento e teste
      data <- as.data.frame(data)
      
      colnames(data) <- c(genotype, env_names, "S2d", "B1", "B1_B2")
      
      data <-
        as_tibble(data) |>
        mutate_at(vars(!!sym(genotype)), as_factor) # colocar para sair com o nome dos ambientes
      
      return(list(
        data1 = data,
        data2 = data[,-c((n_environ + 3):(n_environ + 4))],
        data3 = data[,-c((n_environ + 2):(n_environ + 4))],
        time = time
      ))
      
    } else {
      print("One or more of the following errors has occured:")
      print("The length of the list used in parameters is not 7;")
      print("One or more arguments is not numeric;")
      print("Alpha is too small or too large.")
    }
    
  }


saveRDS(train.set_ctv, "Functions/train.set_ctv.rds")


# Function for graphical plotting ----------------------------------------------------------

#' @param data input (data containing the means of each genotype in each environment);
#' @param jump jump between the visualization of the observations (genotypes);
#' @param y_lab text in y-axis.
#' @param x_lab text in x-axis;
#' @param env_index environmental index. If empty, data will be used to estimate it;
#' @param obs_name name of the observations (example: genotype);
#' @param group_name group that belongs the observation;
#' @param obs_column column in data frame that identifies the group in which each genotype belongs to;
#' @param y_lim limits of the y-axis;
#' @param y_break y-axis scale intervals;
#' @param sec wait in seconds for graphical plots.


show_graphic <-
  function(data,
           jump = 1,
           y_lab = "Production (Kg/ha)",
           env_index = c(),
           x_lab = "Environmental index",
           obs_name = "Genotype",
           group_name = "Group",
           obs_column = 1,
           y_lim = c(2000, 7000),
           y_break = 1000,
           sec = 3)
  {
    g = c()
    if (is.factor(data[, obs_column][[1]])) {
      g <- as.matrix(data[, obs_column])
      data <- as.matrix(data[, -obs_column])
    } else {
      data <- as.matrix(data)
    }
    
    if (length(env_index) != 0) {
      env_index = env_index
    } else {
      env_index = colMeans(data) - mean(as.matrix(data))
    }
    
    
    
    for (i in seq(1, dim(data)[1], jump)) {
      grafico = ggformula::gf_point(data[i, ] ~ env_index,
                                    ylab = y_lab,
                                    xlab = x_lab,
                                    alpha = 0.5) +
        ggplot2::ggtitle(paste(obs_name, i),
                         subtitle = paste(group_name, g[i])) +
        ggplot2::theme_bw() +
        ggplot2::scale_y_continuous(limits = y_lim,
                                    breaks = seq(y_lim[1], y_lim[2], y_break))
      
      print(grafico)
      Sys.sleep(sec)
      if (i == dim(data)[1] - jump + 1)
        break
      
    }
  }

saveRDS(show_graphic, "Functions/show_graphic.rds")


# Function to train the neural network ----------------------------------------------------

#' Par?metros da fun??o:
#' @param data dados de treinamento e valida??o.
#' @param decay taxa de decaimento.
#' @param maxit m?ximo de itera??es do treinamento.
#' @param erro.max erro desejado.
#' @param neurons sequ?ncia de onde ser?o selecionados os neur?nios das camadas ocultas.
#' @param layers n?mero de camadas ocultas.
#' @param column_group coluna que corresponde aos grupos
#' @param p_train valor que determina a quantidade dos dados utilizados para treinamento


train.mlp.net <-
  function(data,
           decay = 0.0005,
           maxit = 5000,
           error.max = 0.05,
           neurons = seq(10, 30),
           layers = 1,
           column_group = 1,
           p_train = 0.8,
           norm_type = "norm")
  {
    len = max(as.numeric(levels(data[, column_group][[1]])))
    
    n_grupos = dim(data)[1] / len # encontra o n?mero de observa??es em cada grupo
    n_t = n_grupos * p_train # no. de obs que ser?o utilizadas para treinamento em cada grupo
    i <- 1
    size <- matrix(NA, 10000, layers)
    erro <- matrix(NA, 10000)
    repeat {
      ad = 0
      samp = c()
      for (z in 1:len) {
        set <- sample(1:n_grupos, n_t, replace = FALSE) + ad * n_grupos
        ad = ad + 1
        samp[(n_t * (z - 1) + 1):(n_t * (z))] = set # linhas a semre selecionadas para o treinamento
      }
      
      dadosX <- as.matrix(data[, -column_group])
      dadosY <- decodeClassLabels(as.matrix(data[, column_group]))
      
      
      # Algoritmo para consertar dadosY
      dadosY2 <-
        matrix(0, length(as.numeric(data[, column_group][[1]])), len)
      
      
      for (k in 1:dim(dadosY)[2]) {
        for (l in 1:dim(dadosY2)[2]) {
          if (as.numeric(colnames(dadosY))[k] == l) {
            dadosY2[, l] = dadosY[, k]
          }
        }
      }
      colnames(dadosY2) <- paste0(1:len)
      
      inputsTrain <- dadosX[samp,]
      targetsTrain <- dadosY2[samp,]
      inputsTest <- dadosX[-samp,]
      targetsTest <- dadosY2[-samp,]
      
      dados1 <-
        list(
          inputsTrain = inputsTrain,
          targetsTrain = targetsTrain,
          inputsTest = inputsTest,
          targetsTest = targetsTest
        )
      
      dados_norm <- normTrainingAndTestSet(dados1, type = norm_type)
      
      
      ncam <- sample(neurons, layers, replace = TRUE)
      size[i, 1:layers] <- ncam
      
      rede_mlp <- mlp(
        x = dados_norm$inputsTrain,
        y = dados_norm$targetsTrain,
        learnFuncParams = decay,
        size = ncam,
        maxit = maxit,
        linOut = FALSE
      )
      
      pred_mlp <- predict(rede_mlp, dados_norm$inputsTest)
      b <- RSNNS::confusionMatrix(dados_norm$targetsTest, pred_mlp)
      
      # Algoritmo para consertar a tabela de confus?o.
      c <- matrix(0, len, len)
      
      
      for (k in 1:dim(b)[2]) {
        for (l in 1:dim(c)[2]) {
          if (as.numeric(colnames(b))[k] == l) {
            c[, l] = b[, k]
          }
        }
      }
      
      dimnames(c) <- list(1:len, 1:len)
      
      # Erro pela tabela de confus?o (c).
      erro[i] = eb = 1 - (sum(diag(c)) / sum(c))
      
      
      # Erro pela tabela de confus?o (c).
      erro[i] = eb = 1 - (sum(diag(c)) / sum(c))
      
      print(ncam)
      print(eb)
      print(c)
      
      if (i == 1) {
        rede.menor.erro <- rede_mlp
        menor.erro <- eb
        matriz_c <- c
        neu.util <- ncam
        metrics <- caret::confusionMatrix(matriz_c)
        indexes <-
          c(metrics$overall[1:2], colMeans(metrics$byClass)[c(1:2, 5:7)])
      }
      if (eb <= menor.erro) {
        rede.menor.erro <- rede_mlp
        menor.erro <- eb
        matriz_c <- c
        neu.util <- ncam
        metrics <- caret::confusionMatrix(matriz_c)
        indexes <-
          c(metrics$overall[1:2], colMeans(metrics$byClass)[c(1:2, 5:7)])
      }
      if (erro[i] < error.max)
        break
      i <- i + 1
    }
    return(
      list(
        mlp_with_smallest_error = rede.menor.erro,
        smallest_error = menor.erro,
        confusion_matrix = matriz_c,
        used_neurons = neu.util,
        all_neurons = size,
        all_erros = erro,
        indexes = indexes
      )
    )
  }


saveRDS(train.mlp.net, "Functions/train.mlp.net.rds")


# Function for classification of adaptability --------------------------------------------

#' @param data organized real data that will be submitted to classification;
#' @param trained.net trained network;
#' @param norm_type type of normalization. Implemented are "0_1", "center", and "norm".


net.adaptability.analysis <-
  function(data, trained.net, norm_type = "norm") {
    g <- as.vector(data_for_classification[, 1][[1]])
    data <- data[,-1]
    last.col <- dim(data)[2]
    data <- as.matrix(data)
    ng <- nrow(data[, -last.col]) # number of genotypes
    na <- ncol(data[, -last.col]) # number of environments.
    
    Both <- round(rowMeans(data[, -last.col]), 2)
    ia <- colMeans(data[, -last.col]) - mean(data[, -last.col])
    Unfavorable <- round(rowMeans(data[, -last.col][, ia < 0]), 2)
    Favorable <- round(rowMeans(data[, -last.col][, ia >= 0]), 2)
    
    data <- normalizeData(data, type = norm_type)
    
    # Application of the trained net to the set of real data
    pred <- predict(trained.net, data, type = "class")
    c <- as.numeric(apply(pred, 1, which.max))
    
    Adaptability <- NULL
    Stability <- NULL
    
    for (i in 1:ng)
    {
      if (c[i] == 1 |
          c[i] == 10)
        Adaptability[i] <- ("(<1 | <1) Exclusive for unfavorable env.")
      if (c[i] == 2 |
          c[i] == 11)
        Adaptability[i] <- ("(<1 | =1) Exclusive for unfavorable env.")
      if (c[i] == 3 | c[i] == 12)
        Adaptability[i] <- ("(<1 | >1) Ideal")
      
      if (c[i] == 4 |
          c[i] == 13)
        Adaptability[i] <- ("(=1 | <1) Not recomended")
      if (c[i] == 5 |
          c[i] == 14)
        Adaptability[i] <- ("(=1 | =1) General")
      if (c[i] == 6 |
          c[i] == 15)
        Adaptability[i] <- ("(=1 | >1) Exclusive for favorable env.")
      
      if (c[i] == 7 |
          c[i] == 16)
        Adaptability[i] <- ("(>1 | =1) Not recomended")
      if (c[i] == 8 |
          c[i] == 17)
        Adaptability[i] <- ("(>1 | <1) Not recomended")
      if (c[i] == 9 |
          c[i] == 18)
        Adaptability[i] <- ("(>1 | >1) Exclusive for favorable env.")
      
      if (c[i] >= 1 & c[i] <= 9)
        Stability[i] <- ("High")
      if (c[i] >= 10 & c[i] <= 18)
        Stability[i] <- ("Low")
      
    }
    means <- data.frame(cbind(Unfavorable, Favorable, Both))
    adap_stab <- data.frame(cbind(data_for_classification[, 1], 
                                  Adaptability, Stability))
    result <- cbind(adap_stab, means)
    
    return(result)
    
  }

saveRDS(net.adaptability.analysis,
        "Functions/net.adaptability.analysis.rds")


# Function for classification visualization 1 ---------------------------------------------

gt_table <- function(ranking) {
  table_gt = ranking %>%
    as.data.frame() %>%
    gt() %>%
    tab_footnote(footnote = "Values greater than the column mean highlighted in red.",
                 locations = cells_column_labels(columns = 4:6)) %>%
    tab_footnote(footnote = "Adaptability based on the method of Cruz, Torres & Vencovsky (1989).",
                 locations = cells_column_labels(columns = 2)) %>%
    tab_footnote(footnote = "Stability based on the method of Finlay & Wilkson (1963).",
                 locations = cells_column_labels(columns = 3)) %>%
    tab_spanner(label = "Average of Genotypes",
                columns = 4:6) %>%
    gt_highlight_rows(
      rows = which(
        ranking$Adaptability == "(<1 | >1) Ideal" &
          ranking$Stability == "High" &
          ranking$Both > mean(ranking$Both)
      ),
      font_weight = "bold"
    ) %>%
    tab_style(
      style = cell_text(weight = "bold", color = "#C42754"),
      locations = cells_body(
        columns = 4,
        rows = Unfavorable >= mean(ranking$Unfavorable)
      )
    ) %>%
    tab_style(
      style = cell_text(weight = "bold", color = "#C42754"),
      locations = cells_body(
        columns = 5,
        rows = Favorable >= mean(ranking$Favorable)
      )
    ) %>%
    tab_style(
      style = cell_text(weight = "bold", color = "#C42754"),
      locations = cells_body(columns = 6,
                             rows = Both >= mean(ranking$Both))
    ) %>%
    cols_align(align = "center",
               columns = c(1, 3:6))
  return(table_gt)
  
}

saveRDS(gt_table, "Functions/gt_table.rds")


# Function for classification visualization 2 ---------------------------------------------

ft_table <- function(ranking) {
  set_flextable_defaults(
    font.size = 10,
    theme_fun = theme_vanilla,
    padding = 6,
    background.color = "white",
    border.color = 'black'
  )
  
  dat_mean <- summarize(ranking,
                        across(Unfavorable:Both, ~ mean(.x, na.rm = TRUE))) %>%
    add_column(Gen = "Mean", .before = TRUE)
  dat <- bind_rows(ranking, dat_mean)
  
  
  ft <- flextable(dat, col_keys = colnames(ranking)) %>%
    colformat_double(digits = 2, j = colnames(ranking)) %>%
    bold(i = nrow_part(., part = "body"))  %>%
    color(
      i = 1:dim(ranking)[1],
      j = 4:6,
      color = function(x)
        ifelse(x > mean(x), "#C42754", "transparent")
    ) %>%
    bold(.,
         ~ Unfavorable > mean(ranking$Unfavorable),
         ~ Unfavorable,
         bold = TRUE) %>%
    bold(., ~ Both > mean(ranking$Both),
         ~ Both, bold = TRUE) %>%
    bold(., ~ Favorable > mean(ranking$Favorable),
         ~ Favorable, bold = TRUE)  %>%
    footnote(
      i = 1,
      j = 2,
      value = as_paragraph(
        c(
          "Adaptability based on the method of Cruz, Torres & Vencovsky (1989)"
        )
      ),
      ref_symbols = "1",
      part = "header",
      inline = TRUE
    ) %>%
    footnote(
      i = 1,
      j = 3,
      value = as_paragraph(c(
        "Stability based on the method of Finlay & Wilkson (1963)"
      )),
      ref_symbols = "2",
      part = "header",
      inline = TRUE
    ) %>%
    footnote(
      i = 1,
      j = 4:6,
      value = as_paragraph(c(
        "Values greater than the average of the column in red."
      )),
      ref_symbols = "3",
      part = "header",
      inline = TRUE
    ) %>%
    padding(j = 2,
            padding.right = 110,
            part = "header") %>%
    add_header_row(
      values = c("", "Average of Genotypes"),
      colwidths = c(3, 3),
      top = T
    ) %>%
    theme_booktabs() %>%
    align(j = c(3:6),
          align = "center",
          part = "all") %>%
    bg(
      i = which(
        ranking$Adaptability == "(<1 | >1) Ideal" &
          ranking$Stability == "High" &
          ranking$Both > mean(ranking$Both)
      ),
      bg = "#FEF79C"
    ) %>%
    bold(i = 1:2, bold = TRUE, part = "header") %>%
    hline(.,
          i = ~ before(Gen, "Mean"),
          border = fp_border_default(width = 2)) %>%
    autofit()
  return(table_ft = ft)
  
}

saveRDS(ft_table, "Functions/ft_table.rds")
