# Mathematical specification

This file is the only mathematical specification for the condition-comparable Pando workflow.

## Candidate dictionary

For one broad cell type, target gene \(g\), pooled cells, and conditions \(c\), let \(E_g^{pool}\) be the pooled candidate set and \(E_{g,c}\) the candidate set found within condition \(c\). The frozen dictionary is the union of complete TF–peak–target triples:

\[
E_g^{\cup}=E_g^{pool}\cup\bigcup_c E_{g,c}.
\]

TFs, peaks, and targets are not unioned independently and are never recombined by Cartesian product.

## Condition fit

For edge \(e\), paired cell \(i\), and condition \(c\), the unscaled predictor is

\[
x_{e,i}=RNA_{TF(e),i}\,ATAC_{peak(e),i}.
\]

Every condition uses the same ordered dictionary and fits

\[
y_{g,i}=\alpha_{g,c}+\sum_{e\in E_g^{\cup}}\beta_{e,g,c}x_{e,i}+\varepsilon_{g,i},
\qquad i\in c.
\]

The implementation uses a Gaussian identity GLM with an intercept, `interaction_term = ":"`, and `scale = FALSE`. Pooled coefficients are not used to centre, scale, or calibrate condition coefficients.

## Estimability and edge selection

Let \(a_{e,g,c}\) indicate that the coefficient is estimable. Zero-variance, aliased, non-finite, and insufficient-residual-degree-of-freedom coefficients remain unavailable.

For the configured adjusted-P threshold \(q\),

\[
s_{e,g,c}=\mathbf{1}\{a_{e,g,c}=1\ \land\ padj_{e,g,c}<q\}.
\]

The effect exported for downstream projection is

\[
\theta_{e,g,c}=
\begin{cases}
\widehat\beta_{e,g,c}, & s_{e,g,c}=1,\\
0, & s_{e,g,c}=0.
\end{cases}
\]

The complete coefficient table is retained. An unavailable coefficient remains `NA`; a non-significant coefficient is not interpreted as a biological zero.

## Paired-cell projection

For cell \(i\) in condition \(c\), the target-level regulatory score is

\[
G_{i,g,c}=\sum_{e\in E_g^{\cup}}\theta_{e,g,c}x_{e,i}.
\]

For aggregation group \(u\) with membership set \(M_u\),

\[
G_{u,g}=\frac{1}{|M_u|}\sum_{i\in M_u}G_{i,g,c(i)}.
\]

Projection is performed on paired cells before aggregation. The fitted coefficient and the RNA and ATAC values always come from the cell's own condition and the recorded preprocessing reference.

## Inference scope

The reported Gaussian GLM P values are conditional on the frozen candidate dictionary. They do not include selective-inference correction for candidate discovery.
