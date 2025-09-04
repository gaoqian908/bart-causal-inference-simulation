# bart-causal-inference-simulation

This repository includes R code for running simulations to compare different causal inference methods, mainly focusing on BART (Bayesian Additive Regression Trees).

## What’s inside

There are two sets of code:

1. One compares four methods:  
   - BART  
   - Outcome Regression  
   - Inverse Probability Weighting (IPW)  
   - Augmented IPW (AIPW)

2. The other compares all of the above plus a BART-based AIPW method.

## About the simulation

- The current code is set to run with `iterations = 1000`.
- If you want to test with `iterations = 100` or `iterations = 500`, just manually change the value of `n` in the script.
- The simulations are based on three different response surfaces with varying complexity and overlap.

## Notes

- I referred to parts of the code from Jennifer Hill’s paper on BART for causal inference.
- This is a personal project done as part of my dissertation, and the code is meant for learning and academic use.

Thanks for checking it out!
