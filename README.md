# Optimizing the dispatch of energy generators 
## Summary
This project utilizes convex optimization with application to the economic dispatch problem in power systems. The purpose of this is to optimally schedule distributed energy generators in a distribution feeder to minimize economic costs, while maintaining safe operating constraints. I utilize the convex DistFlow equations of Baran & Wu, which model power flow through radial distribution power networks. These equations model active & reactive power, branch power flow limits, node voltage limits, etc. This was completed as part of my Energy Systems Control class in my MS program. Thanks to Prof. Scott Moura for providing data, overall framework and cvxpy syntax support for this problem. 

The mathematical formulation of this problem is not included in this document, although it is necessary to understand the constraints and forumlas presented in this notebook.

## Key questions explored
- Can we dispatch generators in a electrical network to minimize cost?
- How does the solution change when considering variable renewable generation?

## Techniques used
- convex optimization and relaxation
- energy systems engineering
- second-order cone programming to deal with stochasticity

## Key findings
- As we increase the complexity of the model, we add more constraints to power flow. This increases the minimum cost.
- The cost when considering renewable generation is slightly than if there were no renewables. This makes sense because renewables are variable and we don't know exactly how much power will be available from them. Thus, if we include renewables in our power flow schedule, we’ve introduced uncertainty and thus risk. There is a cost associated with this risk.
