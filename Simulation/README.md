This folder includes the materials for conducting Monte Carlo Simulations to evaluate our proposed Huber loss-based super learber.

We evaluated the performance of the proposed method using Monte Carlo simulations. 
We were interested in assessing performance of the method in two tasks that are common in biomedical cost research: cost prediction and 
estimation of effects of interventions on costs. In both cases, we focus on the use case of our method where the Huber loss is used as a 
\emph{working} loss function, as we anticipate these situations may arise more readily in practice. 
In both simulation studies, we used the same data generating process for $X$. 
The vector $X$ consisted of ten components, distributed as $X_1 \sim \mbox{Bernoulli}(0.5)$, 
$X_2 \sim \mbox{Uniform}(0,1)$, $X_3 \sim \mbox{Normal}(0,1)$, $X_4 \sim \mbox{Gamma}(1,1)$, $X_5 \sim \mbox{Poisson}(1)$, 
$X_6 \sim \mbox{Bernoulli}(0.2)$, $X_7 \sim \mbox{Uniform}(-1,1)$, $X_8 \sim \mbox{Normal}(0,3)$, 
$X_9 \sim \mbox{Gamma}(0.5,1)$, $X_{10} \sim \mbox{Poisson}(2)$. In both simulations, the variables $X_1, \dots, 
X_5$ impacted the distribution of costs, while the others were noise. The other aspects of the data generating process differed by 
simulation setting and are described below.
