from pathlib import Path


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected one match, found {count}")
    return text.replace(old, new, 1)

# Remove the obsolete condition_grn_fit(network_name=) compatibility surface.
path = Path("R/condition_grn.R")
text = path.read_text()
text = replace_once(
    text,
    "#' @param network_name Retained for compatibility; fit selection is by cell type.\n",
    "",
    "condition_grn_fit roxygen compatibility parameter",
)
old = '''condition_grn_fit.GRNData <- function(
    object, network_name = NULL, cell_type = NULL, ...) {
    fits <- object@grn@params$condition_grn_fits
'''
new = '''condition_grn_fit.GRNData <- function(
    object, cell_type = NULL, ...) {
    dots <- list(...)
    if (length(dots)) {
        label <- names(dots)
        label[!nzchar(label)] <- "<unnamed>"
        stop(
            "Unused condition_grn_fit argument(s): ",
            paste(label, collapse = ", "), call. = FALSE
        )
    }
    fits <- object@grn@params$condition_grn_fits
'''
text = replace_once(text, old, new, "condition_grn_fit signature")
path.write_text(text)

path = Path("man/condition_grn_fit.Rd")
text = path.read_text()
text = replace_once(
    text,
    "\\method{condition_grn_fit}{GRNData}(object, network_name = NULL, cell_type = NULL, ...)\n",
    "\\method{condition_grn_fit}{GRNData}(object, cell_type = NULL, ...)\n",
    "condition_grn_fit Rd usage",
)
text = replace_once(
    text,
    "\\item{network_name}{Retained for compatibility.}\n",
    "",
    "condition_grn_fit Rd parameter",
)
text = text.replace(
    "\\item{...}{Additional arguments.}",
    "\\item{...}{Must be empty. Unknown arguments are rejected.}",
)
path.write_text(text)

# Retire the obsolete shared-baseline tutorial, which references removed APIs.
obsolete = Path("docs/condition-subgrn-shared-baseline.md")
if obsolete.exists():
    obsolete.unlink()

# Make the current release contract explicit without rewriting historical NEWS.
path = Path("NEWS.md")
text = path.read_text()
heading = "# Pando 2.0.0\n"
if not text.startswith(heading):
    current = '''# Pando 2.0.0

- Replaces the retired condition-aware nested-CV, sparse-group, native solver,
  structural-zero and OOF projection engines with one two-stage
  common-dictionary workflow.
- Discovers exact TF-peak-target candidates in the complete cell type and in
  each condition, freezes their exact union, and fits the same Gaussian
  identity `TF:peak` design in every condition with `scale = FALSE`.
- Records RNA/ATAC layers, ATAC value semantics and a preprocessing fingerprint
  in candidate tables, the frozen dictionary, fitted networks and projections.
  Candidate unions and fixed fits reject mixed preprocessing references.
- Uses BH-adjusted condition coefficients for `penalty_effect`; unavailable or
  aliased coefficients remain unavailable rather than becoming fitted zeros.
- Removes obsolete condition APIs and compatibility arguments. The supported
  condition API is `discover_grn_edges()`, `union_grn_edges()`,
  `fit_grn_from_edges()`, `infer_condition_grn()`, `condition_grn_fit()`,
  `condition_grn_subgraph()`, `project_condition_grn_cells()` and
  `aggregate_condition_grn_projection()`.

'''
    path.write_text(current + text)

# Focused regression assertions for the supported API only.
path = Path("tests/testthat/test-common-dictionary-direct-source.R")
text = path.read_text()
addition = r'''

test_that("condition fit extraction has no compatibility arguments", {
  method <- getS3method("condition_grn_fit", "GRNData")
  expect_false("network_name" %in% names(formals(method)))
  expect_error(
    method(methods::new("GRNData"), network_name = "legacy"),
    "Unused condition_grn_fit argument"
  )
})

test_that("retired condition APIs are absent", {
  namespace <- asNamespace("Pando")
  retired <- c(
    "condition_grn_contrast",
    "project_condition_grn_primary_cells",
    "project_condition_grn_to_cells"
  )
  expect_false(any(vapply(retired, exists, logical(1), envir = namespace,
                          inherits = FALSE)))
})
'''
if "condition fit extraction has no compatibility arguments" not in text:
    path.write_text(text.rstrip() + addition + "\n")

# Static repository contract.
source = Path("R/condition_grn.R").read_text()
if "network_name = NULL, cell_type = NULL" in source:
    raise RuntimeError("obsolete condition_grn_fit network_name remains")
if Path("docs/condition-subgrn-shared-baseline.md").exists():
    raise RuntimeError("obsolete shared-baseline tutorial remains")
if not Path("NEWS.md").read_text().startswith("# Pando 2.0.0"):
    raise RuntimeError("current NEWS contract is missing")
