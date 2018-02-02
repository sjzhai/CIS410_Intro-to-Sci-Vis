Our goal is to advance 1s.
RK4: h = 0.01
100 steps
4 velocity evaluations per step 
--> 400 velocity evaluations 

RK4 with h of 0.01 --> Euler with h of 1e-8
To adcanve 1s with a step size of 1e-8.
To adcanve 1s with a step size of 0.1. --> 10 (1/0.1)
To adcanve 1s with a step size of 1e-1. --> 10 (1/0.1)
To adcanve 1s with a step size of 1e-8. --> 1/1e-8 = 100M steps
*AND*
1 velocity evaluation perstep
So
100M steps * 1 velocity evaluation per step = 100M velocity evaluations

What is the ratio of velocity evaluations between Euler and RK4 to achieve the same accuracy?
(with an h of 0.01)
-->
Assume advance 1s.
(Assume advancing 0.01s or 100s.)
