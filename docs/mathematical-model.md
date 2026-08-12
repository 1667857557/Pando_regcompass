# Mathematical specification

This file is the mathematical specification for the condition-comparable Pando workflow.

## 1. Structural candidate space

For target gene \(g\), Pando first restricts candidate regulatory peaks by the supplied regulatory regions, the target regulatory domain, and TF motif support. In the RegCompass workflow the supplied region prior is normally

\[
R_{prior}=R_{phastCons}\cup R_{SCREEN\ cCRE}.
\]

For an exact edge \(e=(g,f,r)\), candidate identity always means the complete target-TF-region triple. TFs, peaks, and targets are never unioned independently and recombined by Cartesian product.

## 2. Pando correlation support and frozen common dictionary

For condition \(c\), define

\[
\rho^{peak}_{e,c}=cor_c(ATAC_r,RNA_g),\qquad
\rho^{TF}_{e,c}=cor_c(RNA_f,RNA_g).
\]

Pando uses strict threshold comparisons in the implementation. An edge has local Pando support in condition \(c\) when both configured gates are satisfied in that same condition:

\[
L_{e,c}=\mathbf 1\left(|\rho^{peak}_{e,c}|>\tau_{peak}\right)
\mathbf 1\left(|\rho^{TF}_{e,c}|>\tau_{TF}\right).
\]

The same correlations are also calculated on all eligible-condition cells pooled within the cell type. Let

\[
G_e=\mathbf 1\left(|\rho^{peak}_{e,global}|>\tau_{peak}\right)
\mathbf 1\left(|\rho^{TF}_{e,global}|>\tau_{TF}\right).
\]

The defaults are \(\tau_{TF}=0.05\) and \(\tau_{peak}=0.05\), but these are adjustable parameters rather than fixed parts of the statistical model. The same user-supplied values are used for pooled/global and condition-specific candidate discovery and are stored with the fit.

The frozen common dictionary is

\[
D_g=\operatorname{unique}\left(
D_{g,global}\cup\bigcup_c D_{g,c}
\right)
=\{e:G_e=1\ \lor\ \exists c\ L_{e,c}=1\}.
\]

Deduplication is performed on the exact \((target,TF,region)\) key. Global support is included to reduce false-negative candidate exclusion caused by unstable correlations in small conditions; condition-local support is retained separately as provenance.

The correlation gates have one role only: **dictionary admission**. Once an exact edge belongs to \(D_g\), it is retained as a predictor in every fitted condition regardless of whether that condition itself passed the marginal-correlation gates. Therefore a condition with a noisy marginal correlation can still recover a supported conditional ridge effect for an edge admitted globally or by another condition.

Because the same RNA/ATAC observations are used for correlation screening and coefficient estimation, subsequent P values are conditional on a data-selected dictionary and are not formal selective-inference P values.

## 3. Pando TF-by-ATAC predictor

For edge \(e\), paired cell \(i\), and condition \(c\), the predictor is the original Pando regulatory interaction

\[
x_{eic}=RNA_{TF(e),i}\,ATAC_{peak(e),i}.
\]

For target \(g\), every condition uses the same ordered columns \(e\in D_g\):

\[
y_{g,c}=\alpha_{g,c}\mathbf 1+X_{g,c}\beta_{g,c}+\varepsilon_{g,c}.
\]

This common coordinate system is what makes \(\beta_{g,c,e}\) directly comparable across conditions.

## 4. Common scaling and no-fusion ridge

The model contains a separate intercept for every condition. Between-condition shifts in mean \(TF\times ATAC\) therefore do not identify a within-condition slope. For predictor \(e\), the common scale is the equal-condition RMS within-condition standard deviation:

\[
s_e^2=\frac{1}{K}\sum_{c=1}^{K}
\frac{1}{n_c}\sum_{i\in c}(x_{eic}-\bar x_{ec})^2.
\]

Cross-validation relearns this scale inside each training fold. Exported coefficients are back-transformed to the original TF-by-ATAC units.

For one target and \(K\) conditions, Pando minimizes

\[
\frac{1}{K}\sum_{c=1}^{K}
\frac{\|y_c-X_c\beta_c\|_2^2}{n_c}
+\lambda\sum_{c=1}^{K}\|\beta_c\|_2^2.
\]

There is **no cross-condition fusion penalty**. Algebraically the penalized normal system is block diagonal apart from the shared preprocessing/tuning conventions:

\[
\left(X_c^TX_c+\lambda I\right)^{-1}X_c^Ty_c
\]

up to the equal-condition weighting and the implementation's common multiplicative normalization. Ridge keeps the system identifiable under severe collinearity or raw rank deficiency, but does not make individual coefficients uniquely attributable when several highly correlated TF-peak predictors encode nearly the same signal.

One target-specific \(\lambda_g\) is selected by condition-stratified cross-validation and used for every condition of that target. The default one-standard-error rule chooses the largest lambda within one standard error of minimum equal-condition validation MSE.

## 5. Approximate ridge inference and active condition GRNs

Let \(a_{e,c}\) denote that edge \(e\) is estimable in condition \(c\). Robust sandwich covariance for the penalized estimator is used to form approximate ridge-Wald statistics. P values are BH-adjusted separately within each condition.

Define condition-specific statistical support

\[
S_{e,c}=\mathbf 1\{a_{e,c}=1\land q_{e,c}<\alpha\}.
\]

Because every fitted coefficient row already corresponds to an edge in the frozen common dictionary, the active condition GRN is simply

\[
\boxed{A_{e,c}=S_{e,c}}.
\]

The pooled/global and condition-local correlation flags \(G_e\) and \(L_{e,c}\) are retained as provenance describing **why the edge entered the common dictionary**. They are not applied again as a post-fit veto. In particular, if edge \(e\) entered through pooled/global support or through condition \(A\), condition \(B\) may still mark the edge active when \(\widehat\beta_{e,B}\) is estimable and BH-supported, even when \(G_e=0\) and \(L_{e,B}=0\). This is intentional: the marginal-correlation screen and the joint ridge coefficient are different statistics, and reapplying the marginal screen after fitting would reintroduce the small-condition false-negative problem that the common dictionary is meant to avoid.

The downstream regulatory effect is

\[
penalty\_effect_{e,c}=A_{e,c}\,\widehat\beta_{e,c}.
\]

The complete \(\widehat\beta\) table is retained even when an edge is inactive. An inactive coefficient is not interpreted as a biological zero.

## 6. Differential coefficients

Because the dictionary, column identities, scaling convention, and target-specific lambda are shared, direct contrasts are defined on the same coefficient coordinate:

\[
\Delta\beta_e=\beta_{e,A}-\beta_{e,B}.
\]

Pairwise contrasts use the covariance of the no-fusion ridge estimator to form approximate Wald statistics and BH-adjusted contrast P values. Significance in one condition and nonsignificance in another is not used as a substitute for the direct contrast.

## 7. Pando paired-cell projection

For paired cell \(i\) in condition \(c\), Pando's cell-level target regulatory score is

\[
G_{i,g,c}=\sum_{e\in D_g}
penalty\_effect_{e,c}\,RNA_{TF(e),i}\,ATAC_{peak(e),i}.
\]

For aggregation group \(u\) with membership set \(M_u\), the Pando projection helper can aggregate these cell-level scores as

\[
G_{u,g}=\frac{1}{|M_u|}\sum_{i\in M_u}G_{i,g,c(i)}.
\]

RegCompass Layer 1 uses a distinct downstream metacell estimand, \(\beta\,\overline{TF}\,\overline{ATAC}\), documented in the RegCompass mathematical specification; it must not be confused with \(\beta\,\overline{TF\times ATAC}\).

## 8. Inference scope

The coefficient and contrast P values are approximate and conditional on:

1. structural region/domain/motif candidate restrictions;
2. pooled/global or condition-wise Pando correlation screening;
3. the resulting frozen exact-edge dictionary;
4. CV-selected ridge lambda.

They are not exact post-selection inference and are not donor-level population inference when cells rather than biological replicates are the statistical units.
