# Introductory code to simulate a COVID-19 outbreak
# Enrique M Saldarriaga
# The CHOICE Institute, University of Washington
# Jul 2020
# Original post: <URL>

# The goal of this R-script is to simulate an COVID-19 outbreak in a hypothetical territory using real-world parameters

### 1. Defining the infectious model
# Here we need to define the "version" of the model we want to use
# Note the following state variables:
## S - number of susceptible in the population
## E - number of infected individuals in the population
## I - number of infectious individuals in the population
## R - number of recovered/immune individuals in the population
## N - the total population size; N = S + E + I + R
## The most basic model includes S and I states, and the most complex would include all of them
## The main difference between E and I is that exposed people are infected but do not transmit the virus
## we differentiate E and I when there's a known latency period between infection and infectiousness
# For this example, we would use a SEIR model. For a given point in time, the state variables would be defined as:
## S(t+1) = S(t) - lambda*S(t)
## E(t+1) = E(t) + lambda*S(t) - sigma*E(t)
## I(t+1) = I(t) + sigma*E(t) - gamma*I(t)
## R(t+1) = R(t) + gamma*I(t)
## N(t+1) = S(t+1) + E(t+1) + I(t+1) + R(t+1)
## The model is governed by lambda, sigma, and gamma. 
## Respectively, the transition probability of becoming infected, infectious, and recovered
## lambda = beta * I/N
## beta = contact_rate * infectivity; infectivity aka rho
## sigma = 1/latency_period
## gamma = 1/infectious_period; aka duration of the disease

### 2. Implement the model: Use of ODE
# We use Ordinary Differential Equations to describe the discrete change in the state variables over time
# Even with very small time-changes (e.g. 0.1*day, 1 hour), the model is still discrete
# However, smaller changes provide "smoother" lines by observing the change more frequently
# To do this, we use the first 4 equations from above and take the partial derivative over time
## dSdt = -beta*S*I/N; read: change in S over time is equal to...
## dEdt = beta*S*I/N - sigma*E
## dIdt = sigma*E - gamma*I
## dRdt = gamma*I

# The following function describes the change in each state variable
seir <- function(t,init,params){
  with(c(as.list(init),params),{
    dSdt = -beta*S*I/N
    dEdt = beta*S*I/N - sigma*E
    dIdt = sigma*E - gamma*I
    list(c(dSdt,dEdt,dIdt))
  })
}
# Note that we didn't define dRdt. This is because we model a closed population (i.e. no deaths, no births)
## so, R is identified by knowing N = S+E+I+R with N constant.

# The function sir uses three parameters:
## t: time change (e.g. days)
## init: initial numbers for the state variables
## params: beta, sigma, gamma

time = 0
N0 = 1000000
propS = 0.999 #proportion of total population susceptible to disease
initial.Pop <- c(S = propS*N0,
                 E = 0, 
                 I = 1 # The outbreak need at least one initial infected to spread the disease
                      ) 
# We decide for I = 1 to simulate the beginning of the outbreak
# Since N=1,000,0000 = 0.999*1,000,000 + 0 + 1 + R, we are implicitly assuming that R = 999, immune people

initial.Params <- c(beta = 1.5,        # Transmission coefficient
                    sigma = 1/4.2,        # Incubation/Latency period (days)
                    gamma = 1/20,       # Recovery rate: 1 / infectious period (days)
                    N = N0)            # Population size (constant)

### 3. Try the function
# Use the function to observe the estimate changes in a single point of time
seir(t=0,init=initial.Pop,params=initial.Params)
# This produces a list with values: -1.4985  1.4985 -0.0500 
# This means that at time 0, 1.5 (1.5*1/1,000,000) came out the S state, the same that went to E
# Similarly, 0.05 (1*1/20) came out of I, and therefore got into R 

### 4. Solve the ODE model
library(deSolve) #install.packages("deSolve")
# We use the [deSolve package](https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf)
## to solve the ODE model that we have created
## Before going there, we need to define the time points of evaluation

mod_time = seq(0,365,0.5) #first, last, step. 
# We want to model a year of the epidemic, retrieving the results every half day

init_SEIR_model=data.frame(ode(
  y = initial.Pop,
  times = mod_time,
  func = seir,    
  parms = initial.Params
))
init_SEIR_model$R = N0 - rowSums(init_SEIR_model[,c("S","E","I")])
# Change format to get rid of scientific notation and decimals
init_SEIR_model[,2:5]=lapply(format(init_SEIR_model[,2:5], digits=c(0), scientific=F), 
                             function(x) as.numeric(as.character(x)))

### 4. Plot the outbreak
# We're going to use ggplot. You can plot each line using geom_line(), but it's much easier to
## reshape the data to long and then plot it
library(reshape2)
init_SEIR_model=melt(init_SEIR_model,id="time")
colnames(init_SEIR_model) = c("Time", "State", "Count")
## Now we plot
library(dplyr); library(ggplot2);library(ggthemes); library(directlabels)

init_SEIR_model %>%
  mutate(infectious = (State == 'I')) %>% #Single out Infectious evolution to observe better
  ggplot(aes(x=Time, y=Count/1E3,color=State))+
  geom_line(aes(linetype=infectious), size=1, alpha=0.8) +
  labs(title = "SEIR Model: COVID-19 Outbreak", 
       subtitle = "Evolution of State Variables over time", 
       x="Day", y="Count (thousand)") +
  theme_bw(base_size = 25) + 
  scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  theme(axis.title = element_text()) +
  scale_linetype_manual(values=c("dashed", "solid"), guide="none") +
  scale_color_manual(values=c('#134678','#C46D0A','#A32123','#137827'))

### 5. Observations over the results of the model

# First, we observe the typical behavior of an outbreak with a peak around day 37.5 at 625,095 infections
init_SEIR_model[init_SEIR_model$Count==max(init_SEIR_model[init_SEIR_model$State=="I",]$Count),]$Time
with(init_SEIR_model[init_SEIR_model$State=="I",],max(Count))
# The reason for this is the rapid depletion of Susceptible given that it's a very infectious disease.
# Note how steep is the shrinkage of Susceptible, alongside with a rapid increase of Exposed, and even faster increase on Infectious
# In fact, people do not spend much time infected because the latency period is only 4.2 days.
# One interesting "invisible" assumption is that Recovered have immunity - then they won't go back to the susceptible compartment.
# As a consequence, once an important proportion of the population is infected (and recovered) there are no more
## susceptible to infect and the outbreak starts to fade

# The basic reproductive number 
# This is number of new infections that one infected can cause in an entirely susceptible population
# It's calculated the product of the likelihood of transmission per contact (aka transmission rate) (rho), contact rate (c) and the duration of the disease
# Recall that we used beta = rho * c. Also, the duration of the disease is the inverse of the recovery rate: 1/gamma
# Hence R0 = 30:
R0= initial.Params["beta"] * 1/initial.Params["gamma"]
# This is a very high number, and another reason why the epidemic progresses so fast.

### 6. Scenario: 'Flattening the curve'

# All social distancing measures had the objective of flattening the curve. But, how?
# Reducing beta. We want to reduce the number of likely infections that one infected would create in a single point of time
## i.e. despite the how long the disease lasts.
# Beta has two components but only the contact rate is susceptible to change due to behavioral changes.
## There's a lot still unknown about the viral load needed to infect someone, and therefore targeting the rho is less realistic.

# Let's model a new outbreak in which we effectively reduce the contact rate (using social distancing, and quarantine)
## and as consequence the beta is 0.6 instead of 1.5. With everything else constant, let's see the new outbreak

new.Params <- c(beta = 0.6,        # Transmission coefficient
                sigma = 1/4.2,        # Incubation/Latency period (days)
                gamma = 1/20,       # Recovery rate: 1 / infectious period (days)
                N = N0)            # Population size (constant)

new_SEIR_model=data.frame(ode(
  y = initial.Pop,
  times = mod_time,
  func = seir,    
  parms = new.Params #only change in the function definition
))
new_SEIR_model$R = N0 - rowSums(new_SEIR_model[,c("S","E","I")])
new_SEIR_model[,2:5]=lapply(format(new_SEIR_model[,2:5], digits=c(0), scientific=F), 
                            function(x) as.numeric(as.character(x)))
new_SEIR_model=melt(new_SEIR_model,id="time")
colnames(new_SEIR_model) = c("Time", "State", "Count")

new_SEIR_model %>%
  mutate(infectious = (State == 'I')) %>% #Single out Infectious evolution to observe better
  ggplot(aes(x=Time, y=Count/1E3,color=State))+
  geom_line(aes(linetype=infectious), size=1, alpha=0.8) +
  labs(title = "SEIR Model: New COVID-19 Outbreak", 
       subtitle = "Evolution of State Variables over time", 
       x="Day", y="Count (thousand)") +
  theme_bw(base_size = 25) + 
  scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  theme(axis.title = element_text()) +
  scale_linetype_manual(values=c("dashed", "solid"), guide="none") +
  scale_color_manual(values=c('#134678','#C46D0A','#A32123','#137827'))

# We observe that the peak occurs later, on day 65.5 instead of 37.5, and the highest count of infected is 550,446
## instead of 625,095
with(new_SEIR_model[new_SEIR_model$State=="I",],max(Count))
new_SEIR_model[new_SEIR_model$Count==max(new_SEIR_model[new_SEIR_model$State=="I",]$Count),]$Time
# The new R0 is 12. Less than half than before.
R0_new= new.Params["beta"] * 1/new.Params["gamma"]
