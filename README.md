# python_differential_equations
research work on the numerical solution of first-order differential equations
# Task
Using the Euler and Runge-Kutta methods of the 4th order, find a solution to an ordinary differential equation describing the following physical problem:

A force depending on time t according to the law F = b(τ - t)exp(-t/τ) began to act on a resting particle of mass m at the moment t = 0, where b is a positive constant, τ is the time during which this force acts. Find the momentum of a particle as a function of time.

Before writing a program, it is necessary to compose the differential equation itself in the form dx/ dt = f(x, t), where x is the desired function, t is an independent variable (for example, time), f is a known (analytically specified) dependence.

# Implementation

The values of all input parameters are read from the file in.txt in format:<br />
```
<b paremeter>
<tau parameter>
<number of splittings>
<step size>
<t_0 starting time>
<p_0 pressure at t_0>
```


The solution outputs to the files out1.txt and out2.txt.
The solution means a list of time points and a list of values of the desired function at these time points.
