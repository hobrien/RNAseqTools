# Differential Expression

- Once expression is normalised between samples, testing for between-group differences is in principle straightforward, and t-tests have been used for this purpose, but two factors complicate the analysis:
    - Over-dispersion of variance
        - usually accounted for using a negative binomial distribution
    - Confounding factors
        - usually accounted for with a linear modeling framework

- DESeq2: Likelihood ratio tests or Wald tests
- EdgeR: Liklihood ratio tests
