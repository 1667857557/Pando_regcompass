# Mathematical specification

This file defines the condition-comparable multi-task ridge workflow.

## 1. Candidate union

For one broad cell type and target gene \(g\), let \(E_g^{pool}\) be the pooled Pando candidate set and \(E_{g,c}\) the candidate set discovered in condition \(c\). Candidate discovery uses the Pando regulatory-domain, motif, TF-target correlation, and peak-target correlation rules.

The preliminary dictionary is the exact union of complete TF-peak-target triples:

\[
D_g^{(0)}=E_g^{pool}\cup\bigcup_c E_{g,c}.
\]

TFs, peaks, and targets are never unioned independently or recombined by Cartesian product.

For edge \(e\), paired cell \(i\), and condition \(c\), the raw predictor is

\[
x_{e,i}=RNA_{TF(e),i}\,ATAC_{peak(e),i}.
\]

All conditions use the same ordered dictionary.

## 2. Joint multi-task ridge

For target \(g\), the condition coefficients are estimated jointly. Each condition has its own intercept and coefficient vector. The loss is condition-balanced and contains ordinary ridge shrinkage plus a fusion penalty between condition coefficient vectors:

\[
\widehat{\boldsymbol\beta}
=
\arg\min_{\boldsymbol\beta}
\left[
L(\boldsymbol\beta)
+
N\lambda\,\boldsymbol\beta^T P_\gamma\boldsymbol\beta
\right],
\]

where \(P_\gamma\) contains the identity ridge penalty and the configured condition-fusion penalty. `condition_ridge_control$fusion_ratio` controls the fusion term and \(\lambda\) is chosen by cross-validation (`1se` by default).

Predictors are internally scaled by the equal-condition RMS of within-condition standard deviations. Exported coefficients are transformed back to raw `RNA * ATAC` units, so the condition coefficients for one edge share the same final scale.

## 3. Preliminary screening inference

The first joint fit is performed on the full candidate dictionary \(D_g^{(0)}\) with approximate ridge-Wald inference enabled.

For a condition-edge coefficient,

\[
z_{e,c}^{screen}
=
\frac{\widehat\beta_{e,c}^{screen}}
{SE(\widehat\beta_{e,c}^{screen})},
\]

with the standard error obtained from the regularized sandwich covariance used by the implementation. The two-sided screening P value is

\[
p_{e,c}^{screen}=2\Phi\left(-\left|z_{e,c}^{screen}\right|\right).
\]

Within each condition, finite estimable screening P values are adjusted once by Benjamini-Hochberg:

\[
q_{e,c}^{screen}=BH_c(p_{e,c}^{screen}).
\]

For the configured `padj_threshold = q`, condition-specific screening support is

\[
S_{e,c}
=
\mathbf 1\left\{
\text{estimable}_{e,c}^{screen}
\land q_{e,c}^{screen}<q
\right\}.
\]

An edge enters the final common dictionary if at least one condition supports it:

\[
D_g^{(1)}
=
\left\{
e\in D_g^{(0)}:\sum_c S_{e,c}>0
\right\}.
\]

The BH-adjusted screening result is therefore a dictionary-selection and condition-support rule. It is the only coefficient-level significance calculation in the condition workflow.

## 4. Final common-dictionary effect refit

After \(D_g^{(1)}\) is fixed, all conditions are jointly refit on exactly that same dictionary using the same multi-task ridge estimator and cross-validation procedure.

The final fit is run with coefficient inference disabled:

\[
\widehat\beta_{e,c}^{final}
=
\operatorname{MultiTaskRidgeRefit}(D_g^{(1)}),
\]

but no final \(SE\), Wald statistic, P value, or BH-adjusted P value is calculated.

This separation is intentional. The preliminary fit supplies statistical screening evidence on the broad candidate model. The final refit supplies comparable quantitative condition effects on the selected shared model. A preliminary screening P value must not be interpreted as the P value of \(\widehat\beta^{final}\).

For downstream projection, the active quantitative edge effect is

\[
\theta_{e,c}
=
\begin{cases}
\widehat\beta_{e,c}^{final}, & S_{e,c}=1\text{ and the final coefficient is estimable},\\
0, & \text{otherwise}.
\end{cases}
\]

Thus all conditions are estimated on the same final predictor space, while condition-specific support remains defined by the preliminary BH screen.

## 5. Condition differences

The final joint refit retains directly comparable coefficient differences,

\[
\Delta_{e,a,b}
=
\widehat\beta_{e,a}^{final}
-
\widehat\beta_{e,b}^{final}.
\]

Because final inference is disabled, the final refit reports these as quantitative contrasts only; it does not attach a second post-selection contrast P value.

## 6. Paired-cell projection

For cell \(i\) in condition \(c\), the target regulatory score is

\[
G_{i,g,c}
=
\sum_{e\in D_g^{(1)}}
\theta_{e,c}\,x_{e,i}.
\]

For aggregation group \(u\) with membership set \(M_u\),

\[
G_{u,g}
=
\frac{1}{|M_u|}
\sum_{i\in M_u}G_{i,g,c(i)}.
\]

Projection is performed on paired cells before aggregation. The fitted coefficient and RNA/ATAC values use the recorded preprocessing reference.

## 7. Inference scope

`screen_pval` and `screen_padj` are approximate ridge-Wald screening quantities conditional on the preliminary candidate dictionary, CV-selected regularization, and fusion structure. Candidate discovery itself is data dependent, so these values are not claimed to provide selective-inference-corrected FDR for the complete candidate-generation process.

The final refit makes no coefficient-level inferential claim. Its `estimate` values are post-screen quantitative effects on a shared dictionary and are used for condition-comparable projection.
