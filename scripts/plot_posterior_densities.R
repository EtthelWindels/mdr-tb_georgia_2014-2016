##-------------------------------------------------------------------------
## Create posterior density plots for two-type birth-death analyses
## 2021-09-16 Etthel Windels
##-------------------------------------------------------------------------



# Load libraries ----------------------------------------------------------

library(ggplot2)


# Function to load log files ----------------------------------------------

loadLog <- function(filename, burninFrac=0.1, subsample=NA) {
  
  df_in <- as.matrix(read.table(filename, header=T))
  
  if (burninFrac>0) {
    n <- dim(df_in)[1]
    df_in <- df_in[-(1:ceiling(burninFrac*n)),]
  }
  
  if (!is.na(subsample)) {
    indices <- unique(round(seq(1, dim(df_in)[1], length.out=subsample)))
    df_in <- df_in[indices,]
  }
  
  return(df_in)
}


# Load log files ----------------------------------------------------------

logfile <- loadLog('path_to_combined_logfile', burninFrac=0)
logfile <- as.data.frame(logfile)
logfile$lambdaS <- logfile$R0_base*logfile$becomeUninfectiousRate1
logfile$lambdaR <- logfile$lambdaS*logfile$lambda_ratio
logfile$R0R <- (logfile$lambdaR/logfile$becomeUninfectiousRate2)
logfile$R0ratioR <- logfile$R0R/logfile$R0_base

log_data <- data.frame(R0=c(logfile$R0_base, logfile$R0R), 
                       lambda=c(logfile$lambdaS,logfile$lambdaR),
                       type=c(rep('DS',length(logfile$R0_base)), rep('MDR',length(logfile$R0R))))
log_data$type <- factor(log_data$type, levels=c('DS','MDR'))


# Define prior distributions ----------------------------------------------

prior_R0 <- rlnorm(100000,0,0.5)
prior_lambda <- rlnorm(100000,0,0.5)
prior_delta <- rlnorm(100000,0,0.1)

prior_data <- data.frame(R0=c(logfile$R0_base, prior_R0), 
                     lambda=c(logfile$lambda_ratio, prior_lambda),
                     type=c(rep("posterior",length(logfile$R0_base)), rep('prior',length(prior_R0))))
prior_data$type <- as.factor(prior_data$type)
prior_data$type <- factor(prior_data$type, levels=c('prior','posterior'))

prior_data_delta <- data.frame(delta=c(logfile$becomeUninfectiousRate1,logfile$becomeUninfectiousRate2,prior_delta),
                               type=c(rep('posterior (susceptible)',length(logfile$becomeUninfectiousRate1)),rep('posterior (MDR)',length(logfile$becomeUninfectiousRate2)), rep('prior',length(prior_delta))))
prior_data_delta$type <- as.factor(prior_data_delta$type)
prior_data_delta$type <- factor(prior_data_delta$type,levels=c('prior','posterior (susceptible)','posterior (MDR)'))


# Plot posterior densities ------------------------------------------------

R0 <- ggplot(log_data, aes(R0,fill=type, colour=type)) +
  geom_density(stat='density',alpha=0.5)+
  labs(title=expression(atop(bold("Reproductive number"))), x=expression(R[e]), y='Posterior density')+ 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.text = element_blank(),
        plot.title = element_text(size=25, hjust=0.5)) +
  xlim(0,2)+
  scale_fill_manual(labels=c('susceptible','MDR'),values=c('lightsteelblue3','moccasin'))+
  scale_colour_manual(labels=c('susceptible','MDR'),values=c('lightsteelblue3','moccasin'))

lambda <- ggplot(log_data, aes(lambda,fill=type, colour=type)) +
  geom_density(stat='density',alpha=0.5)+
  labs(title=expression(atop(bold("Transmission rate"))), x='Transmission rate', y='Posterior density')+ 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.text = element_blank(),
        plot.title = element_text(size=25, hjust=0.5)) +
  xlim(0,2)+
  scale_fill_manual(labels=c('susceptible','MDR'),values=c('lightsteelblue3','moccasin'))+
  scale_colour_manual(labels=c('susceptible','MDR'),values=c('lightsteelblue3','moccasin'))

ggsave("R0.svg", plot=R0, path="output_path",scale=1 )
ggsave("lambda.svg", plot=lambda, path="output_path",scale=1 )



# Plot prior + posterior densities ----------------------------------------

R0prior <- ggplot(prior_data, aes(R0,fill=type, colour=type)) +
  geom_density(stat='density',alpha=0.5)+
  labs( x=expression(R[e]~susceptible), y='Density')+ 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        plot.title = element_text(size=25, hjust=0.5)) +
  xlim(0,3.0)+
  scale_fill_manual(labels=c('prior','posterior'),values=c('lightgrey','darkolivegreen3'))+
  scale_colour_manual(labels=c('prior','posterior'),values=c('lightgrey','darkolivegreen3'))

lambdaprior <- ggplot(prior_data, aes(lambda,fill=type, colour=type)) +
  geom_density(stat='density',alpha=0.5)+
  labs( x="Transmission rate ratio", y='Density')+ 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        plot.title = element_text(size=25, hjust=0.5)) +
  xlim(0,3.0)+
  scale_fill_manual(labels=c('prior','posterior'),values=c('lightgrey','darkolivegreen3'))+
  scale_colour_manual(labels=c('prior','posterior'),values=c('lightgrey','darkolivegreen3'))

deltaprior <- ggplot(prior_data_delta, aes(delta,fill=type, colour=type)) +
  geom_density(stat='density',alpha=0.5)+
  labs( x="Becoming uninfectious rate", y='Density')+ 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        plot.title = element_text(size=25, hjust=0.5)) +
  xlim(0.5,1.5)+
  scale_fill_manual(labels=c('prior','posterior (susceptible)','posterior (MDR)'),values=c('lightgrey','darkolivegreen3','darkolivegreen4'))+
  scale_colour_manual(labels=c('prior','posterior (susceptible)','posterior (MDR)'),values=c('lightgrey','darkolivegreen3','darkolivegreen4'))

ggsave("R0prior.svg", plot=R0prior, path="output_path",scale=1 )
ggsave("lambdaprior.svg", plot=lambdaprior, path="output_path",scale=1 )
ggsave("deltaprior.svg", plot=deltaprior, path="output_path",scale=1 )

