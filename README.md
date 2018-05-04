EE 511 Simulation Methods for Stochastic Systems Project #5: Optimization & Sampling via MCMC
===============================================================================================

[MCMC for Sampling]
-------------------
The random variable X has a mixture distribution: 60% in a Beta(1,8) distribution and 40% in a Beta(9,1)
distribution.

i. Implement a Metropolis-Hastings algorithm to generate samples from this distribution.

ii. Run the algorithm multiple times from different initial points. Plot sample paths for the algorithm.
Can you tell if/when the algorithm converges to its equilibrium distribution?

Plot sample paths for the algorithm using different proposal pdfs. Comment on the effect of lowvariance
vs high-variance proposal pdfs on the behavior of your algorithm.

[MCMC for Optimization]
--------------------------
The n-dimensional Scwefel function
f(x)=418.9829 n−Σ
i=1
n
xisin √|xi|
xi∈[−500,500 ]
is a very bumpy surface with many local critical points and one global minimum. We will explore the
surface for the case n=2 dimensions.

i. Plot a contour plot of the surface for the 2-D surface

ii. Implement a simulated annealing procedure to find the global minimum of this surface

iii. Explore the behavior of the procedure starting from the origin with an exponential, a polynomial, and a
logarithmic cooling schedule. Run the procedure for t={20, 50, 100, 1000} iterations for k=100 runs each.
Plot a histogram of the function minima your procedure converges to.

iv. Choose your best run and overlay your 2-D sample path on the contour plot of the Schwefel function to
visualize the locations your optimization routine explored.

[Optimal Paths]
----------------
The famous Traveling Salesman Problem (TSP) is an NP-hard routing problem. The time to find optimal
solutions to TSPs grows exponentially with the size of the problem (number of cities). A statement of the
TSP goes thus:
“A salesman needs to visit each of N cities exactly once and in any order. Each city is connected to other
cities via an air transportation network. Find a minimum length path on the network that goes through all
N cities exactly once (an optimal Hamiltonian cycle).”

A TSP solution ⃗c=(c1, ... cN ) is just an ordered list of the N cities with minimum path length. We will be
exploring MCMC solutions to small and larger scale versions of the problem.

i. Pick N=10 2-D points in the [0,1000]x[0,1000] rectangle. These 2-D samples will represent the
locations of N=10 cities.
1. Write a function to capture the objective function of the TSP problem:
D(⃗c)=Σ
i =1
N−1
‖ci+1−ci‖

2. Start with a random path through all N cities ⃗c0 (a random permutation of the cities), an initial
high temperature T0, and a cooling schedule Tk=f (T0 , k ).

3. Randomly pick any two cities in your current path. Swap them. Use the difference between the
new and old path length to calculate a Gibbs acceptance probability. Update the path
accordingly.

4. Update your annealing temperature and repeat the previous city swap step. Run the simulated
annealing procedure “to convergence.”

5. Plot the values of your objective function from each step. Plot your final TSP city tour.
ii. Run the Simulated Annealing TSP solver you just developed for N = {40, 400, 1000} cities. Explore
the speed and convergence properties at these different problem sizes. You might want to play with
the cooling schedules.
