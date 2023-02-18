## clecon

A suite of functions, examples and data sets for conducting classical economic analysis, including 
- the computation of the vector of labor values, the vector of prices of production and the (scalar) uniform rates of profit
- quantifying deviation between the vector of _relative_ labor values and the vector of _relative_ prices of production

# The functions

- `ppsraffa1`: Compute values and prices of production for the circulating capital model using the Sraffian approach
- `ppstdint1`: Compute values and prices of production for the circulating capital model using the Standard Interpretation of Marx's labor theory of value
- `ppstdint2`: Compute values and prices of production for the model with capital stock using the Standard Interpretation of Marx's labor theory of value (with uniform wage rates and no unproductive industries)
- `ppnewint1`: Compute values and prices of production for the circulating capital model using the New Interpretation of Marx's labor theory of value
- `ppnewint2`: Compute values and prices of production for the circulating capital model using the New Interpretation of Marx's labor theory of value allowing for wage differential across industries
- `ppnewint3`: Compute values and prices of production for the circulating capital model using the New Interpretation of Marx's labor theory of value allowing for some unproductive industries
- `ppnewint4`: Compute values and prices of production for the circulating capital model using the New Interpretation of Marx's labor theory of value allowing for wage differential across industries and some unproductive industries
- `ppnewint5`: Compute values and prices of production for the model with capital stock using the New Interpretation of Marx's labor theory of value (uniform wage rates and no unproductive industries)
- `nonreg_tests`: Quantify deviation between relative values and relative prices of production using non-regression-based measures
- `reg_tests`: Quantify deviation between relative values and relative prices of production using regression-based measures

The output from each function is a list. 