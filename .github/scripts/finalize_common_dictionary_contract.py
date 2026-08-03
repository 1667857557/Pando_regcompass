from pathlib import Path
import re

ROOT = Path('.')
R_DIR = ROOT / 'R'


def function_span(text, name):
    pattern = re.compile(r'(?m)^' + re.escape(name) + r'\s*<-\s*function\b')
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise RuntimeError(f'expected one definition of {name}, found {len(matches)}')
    start = matches[0].start()
    open_brace = text.find('{', matches[0].end())
    if open_brace < 0:
        raise RuntimeError(f'opening brace not found for {name}')
    depth = 0
    quote = None
    escaped = False
    comment = False
    i = open_brace
    while i < len(text):
        ch = text[i]
        if comment:
            if ch == '\n':
                comment = False
            i += 1
            continue
        if quote is not None:
            if escaped:
                escaped = False
            elif ch == '\\':
                escaped = True
            elif ch == quote:
                quote = None
            i += 1
            continue
        if ch == '#':
            comment = True
        elif ch in ('"', "'"):
            quote = ch
        elif ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                return start, i + 1
        i += 1
    raise RuntimeError(f'closing brace not found for {name}')


def extract_function(text, name):
    start, end = function_span(text, name)
    return text[start:end]


def replace_function(text, name, replacement):
    start, end = function_span(text, name)
    return text[:start] + replacement.rstrip() + text[end:]


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f'{label}: expected one match, found {count}')
    return text.replace(old, new, 1)


condition_path = R_DIR / 'condition_grn.R'
validation_path = R_DIR / 'zz_common_dictionary_validation.R'
projection_path = R_DIR / 'zz_common_dictionary_projection.R'
if not condition_path.exists() or not validation_path.exists() or not projection_path.exists():
    raise RuntimeError('expected common-dictionary source files are missing')

condition = condition_path.read_text()
validation = validation_path.read_text()

# 1. Put the complete biological dictionary validator in the canonical source.
condition = replace_function(
    condition,
    '.condition_validate_dictionary',
    extract_function(validation, '.condition_validate_dictionary')
)

# 2. Apply the residual-df rule directly in the canonical target fitter.
old_status = '''    status <- if (fit$df.residual < as.integer(min_residual_df)) {
        "insufficient_df"
    } else if (any(template$aliased)) {
        "rank_deficient"
    } else {
        "ok"
    }
    list(
'''
new_status = '''    status <- if (fit$df.residual < as.integer(min_residual_df)) {
        "insufficient_df"
    } else if (any(template$aliased)) {
        "rank_deficient"
    } else {
        "ok"
    }
    if (identical(status, "insufficient_df")) {
        template$estimate[] <- NA_real_
        template$std_err[] <- NA_real_
        template$statistic[] <- NA_real_
        template$pval[] <- NA_real_
        template$estimable[] <- FALSE
        template$aliased[] <- TRUE
        coefficient[["(Intercept)"]] <- NA_real_
    }
    list(
'''
condition = replace_once(condition, old_status, new_status, 'residual-df contract')

# 3. Preserve the standard Network fit schema without a late wrapper.
old_fit_bind = '''    coefficient <- do.call(rbind, lapply(result, `[[`, "coefs"))
    fit_table <- do.call(rbind, lapply(result, `[[`, "gof"))
    rownames(coefficient) <- rownames(fit_table) <- NULL
'''
new_fit_bind = '''    coefficient <- do.call(rbind, lapply(result, `[[`, "coefs"))
    fit_table <- do.call(rbind, lapply(result, `[[`, "gof"))
    rownames(coefficient) <- rownames(fit_table) <- NULL
    fit_table$nvariables <- as.integer(fit_table$nvariables_dictionary)
'''
condition = replace_once(condition, old_fit_bind, new_fit_bind, 'Network nvariables alias')

# 4. Record exact input-layer and preprocessing identity in every fit.
helpers = r'''
.condition_hash_object <- function(value) {
    file <- tempfile("pando_condition_fingerprint_", fileext = ".rds")
    on.exit(unlink(file), add = TRUE)
    saveRDS(value, file, version = 2L, compress = FALSE)
    unname(tools::md5sum(file)[[1L]])
}

.condition_preprocessing_fingerprint <- function(
    object, gene_data, peak_data, rna_layer, peak_layer,
    peak_value_type = "normalized") {
    params <- Params(object)
    command_provenance <- tryCatch({
        commands <- object@data@commands
        lapply(commands, function(command) {
            list(
                name = tryCatch(as.character(command@name),
                                error = function(error) NA_character_),
                assay = tryCatch(as.character(command@assay.used),
                                 error = function(error) NA_character_),
                call = tryCatch(as.character(command@call.string),
                                error = function(error) NA_character_)
            )
        })
    }, error = function(error) list())
    .condition_hash_object(list(
        schema = "pando_common_dictionary_preprocessing_v1",
        rna_assay = params$rna_assay,
        peak_assay = params$peak_assay,
        rna_layer = rna_layer,
        peak_layer = peak_layer,
        peak_value_type = peak_value_type,
        gene_matrix_class = class(gene_data),
        peak_matrix_class = class(peak_data),
        gene_dim = dim(gene_data),
        peak_dim = dim(peak_data),
        gene_cells = rownames(gene_data),
        peak_cells = rownames(peak_data),
        gene_features = colnames(gene_data),
        peak_features = colnames(peak_data),
        preprocessing_commands = command_provenance
    ))
}
'''
anchor = '.condition_common_dictionary_schema <- "pando_condition_grn_common_dictionary_v1"\n'
if helpers.strip() not in condition:
    condition = replace_once(condition, anchor, anchor + helpers + '\n', 'fingerprint helpers')

old_prepare_signature = '''.condition_prepare_common_input <- function(
    object, genes = NULL, peak_to_gene_method = c("Signac", "GREAT"),
    upstream = 100000, downstream = 0, extend = 1000000,
    only_tss = FALSE, peak_to_gene_domains = NULL, verbose = TRUE) {
'''
new_prepare_signature = '''.condition_prepare_common_input <- function(
    object, genes = NULL, peak_to_gene_method = c("Signac", "GREAT"),
    upstream = 100000, downstream = 0, extend = 1000000,
    only_tss = FALSE, peak_to_gene_domains = NULL,
    rna_layer = "data", peak_layer = "data",
    peak_value_type = c("normalized", "probability", "other"),
    verbose = TRUE) {
'''
condition = replace_once(condition, old_prepare_signature, new_prepare_signature, 'prepare signature')
condition = replace_once(
    condition,
    '    peak_to_gene_method <- match.arg(peak_to_gene_method)\n',
    '''    peak_to_gene_method <- match.arg(peak_to_gene_method)
    peak_value_type <- match.arg(peak_value_type)
    valid_layer <- function(value) {
        is.character(value) && length(value) == 1L && !is.na(value) &&
            nzchar(trimws(value))
    }
    if (!valid_layer(rna_layer) || !valid_layer(peak_layer)) {
        stop("RNA and ATAC layer names must be non-empty strings.",
             call. = FALSE)
    }
''',
    'layer validation'
)
prepare_start, prepare_end = function_span(condition, '.condition_prepare_common_input')
prepare = condition[prepare_start:prepare_end]
prepare = replace_once(
    prepare,
    'object, assay = params$rna_assay, layer = "data"',
    'object, assay = params$rna_assay, layer = rna_layer',
    'RNA layer use'
)
prepare = replace_once(
    prepare,
    'object, assay = params$peak_assay, layer = "data"',
    'object, assay = params$peak_assay, layer = peak_layer',
    'ATAC layer use'
)
old_common_subset = '''    gene_data <- gene_data[common_cells, , drop = FALSE]
    peak_data_all <- peak_data_all[common_cells, , drop = FALSE]

    features <- intersect(gene_annot$gene_name, genes)
'''
new_common_subset = '''    gene_data <- gene_data[common_cells, , drop = FALSE]
    peak_data_all <- peak_data_all[common_cells, , drop = FALSE]
    if (identical(peak_value_type, "probability")) {
        observed <- if (inherits(peak_data_all, "sparseMatrix")) {
            peak_data_all@x
        } else {
            as.numeric(peak_data_all)
        }
        if (any(!is.finite(observed)) ||
            any(observed < 0 | observed > 1)) {
            stop("Probability-valued ATAC layers must be finite and in [0, 1].",
                 call. = FALSE)
        }
    }
    preprocessing_fingerprint <- .condition_preprocessing_fingerprint(
        object = object, gene_data = gene_data, peak_data = peak_data_all,
        rna_layer = rna_layer, peak_layer = peak_layer,
        peak_value_type = peak_value_type
    )

    features <- intersect(gene_annot$gene_name, genes)
'''
prepare = replace_once(prepare, old_common_subset, new_common_subset, 'input fingerprint creation')
old_return = '''        region_map = region_map,
        params = params,
        peak_to_gene_method = peak_to_gene_method,
'''
new_return = '''        region_map = region_map,
        params = params,
        rna_layer = rna_layer,
        peak_layer = peak_layer,
        peak_value_type = peak_value_type,
        preprocessing_fingerprint = preprocessing_fingerprint,
        peak_to_gene_method = peak_to_gene_method,
'''
prepare = replace_once(prepare, old_return, new_return, 'prepared provenance')
condition = condition[:prepare_start] + prepare + condition[prepare_end:]

# Public candidate discovery accepts and records the common input representation.
discover_start, discover_end = function_span(condition, 'discover_grn_edges')
discover = condition[discover_start:discover_end]
discover = replace_once(
    discover,
    '    peak_to_gene_domains = NULL, parallel = FALSE, verbose = TRUE) {',
    '''    peak_to_gene_domains = NULL, rna_layer = "data", peak_layer = "data",
    peak_value_type = c("normalized", "probability", "other"),
    parallel = FALSE, verbose = TRUE) {''',
    'discover signature'
)
discover = replace_once(
    discover,
    '    source_type <- match.arg(source_type)\n',
    '    source_type <- match.arg(source_type)\n    peak_value_type <- match.arg(peak_value_type)\n',
    'discover peak type'
)
discover = replace_once(
    discover,
    '''        peak_to_gene_domains = peak_to_gene_domains,
        verbose = verbose
''',
    '''        peak_to_gene_domains = peak_to_gene_domains,
        rna_layer = rna_layer, peak_layer = peak_layer,
        peak_value_type = peak_value_type, verbose = verbose
''',
    'discover prepare forwarding'
)
condition = condition[:discover_start] + discover + condition[discover_end:]

# Fixed-dictionary API exposes but strictly constrains the comparable model.
fit_start, fit_end = function_span(condition, 'fit_grn_from_edges')
fit_fun = condition[fit_start:fit_end]
fit_fun = replace_once(
    fit_fun,
    '''    peak_to_gene_domains = NULL, parallel = FALSE, overwrite = FALSE,
    verbose = TRUE) {''',
    '''    peak_to_gene_domains = NULL, method = "glm",
    family = stats::gaussian(link = "identity"), interaction_term = ":",
    scale = FALSE, rna_layer = "data", peak_layer = "data",
    peak_value_type = c("normalized", "probability", "other"),
    parallel = FALSE, overwrite = FALSE, verbose = TRUE) {''',
    'fit signature'
)
fit_fun = replace_once(
    fit_fun,
    '''    if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
''',
    '''    peak_value_type <- match.arg(peak_value_type)
    family_name <- tryCatch(family$family, error = function(error) NA_character_)
    family_link <- tryCatch(family$link, error = function(error) NA_character_)
    if (!identical(method, "glm") || !identical(family_name, "gaussian") ||
        !identical(family_link, "identity") ||
        !identical(interaction_term, ":") || !identical(scale, FALSE)) {
        stop(
            "Comparable fixed-dictionary fitting requires method='glm', ",
            "Gaussian identity, interaction_term=':', and scale=FALSE.",
            call. = FALSE
        )
    }
    if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
''',
    'fixed model validation'
)
fit_fun = replace_once(
    fit_fun,
    '''        only_tss = only_tss, peak_to_gene_domains = peak_to_gene_domains,
        verbose = verbose
''',
    '''        only_tss = only_tss, peak_to_gene_domains = peak_to_gene_domains,
        rna_layer = rna_layer, peak_layer = peak_layer,
        peak_value_type = peak_value_type, verbose = verbose
''',
    'fit prepare forwarding'
)
condition = condition[:fit_start] + fit_fun + condition[fit_end:]

# Add input provenance to each Network object and combined fit contract.
old_network_params = '''            edge_dictionary = edge_dictionary,
            scale = FALSE,
            interaction = ":",
            adjust_method = adjust_method,
'''
new_network_params = '''            edge_dictionary = edge_dictionary,
            scale = FALSE,
            interaction = ":",
            rna_layer = prepared$rna_layer,
            peak_layer = prepared$peak_layer,
            peak_value_type = prepared$peak_value_type,
            preprocessing_fingerprint = prepared$preprocessing_fingerprint,
            adjust_method = adjust_method,
'''
condition = replace_once(condition, old_network_params, new_network_params, 'Network provenance')

# Multi-condition entry point uses the same explicitly recorded representation.
infer_start, infer_end = function_span(condition, 'infer_condition_grn.GRNData')
infer_fun = condition[infer_start:infer_end]
infer_fun = replace_once(
    infer_fun,
    '''    peak_to_gene_domains = NULL, tf_cor = 0.1, peak_cor = 0,
''',
    '''    peak_to_gene_domains = NULL, rna_layer = "data", peak_layer = "data",
    peak_value_type = c("normalized", "probability", "other"),
    tf_cor = 0.1, peak_cor = 0,
''',
    'infer signature'
)
infer_fun = replace_once(
    infer_fun,
    '    rank_action <- match.arg(rank_action)\n',
    '    rank_action <- match.arg(rank_action)\n    peak_value_type <- match.arg(peak_value_type)\n',
    'infer peak type'
)
infer_fun = replace_once(
    infer_fun,
    '''        only_tss = only_tss,
        peak_to_gene_domains = peak_to_gene_domains,
        verbose = verbose
''',
    '''        only_tss = only_tss,
        peak_to_gene_domains = peak_to_gene_domains,
        rna_layer = rna_layer, peak_layer = peak_layer,
        peak_value_type = peak_value_type, verbose = verbose
''',
    'infer prepare forwarding'
)
infer_fun = replace_once(
    infer_fun,
    '''            rna_assay = prepared$params$rna_assay,
            atac_assay = prepared$params$peak_assay
''',
    '''            rna_assay = prepared$params$rna_assay,
            atac_assay = prepared$params$peak_assay,
            rna_layer = prepared$rna_layer,
            peak_layer = prepared$peak_layer,
            peak_value_type = prepared$peak_value_type,
            preprocessing_fingerprint = prepared$preprocessing_fingerprint
''',
    'fit contract provenance'
)
condition = condition[:infer_start] + infer_fun + condition[infer_end:]
condition_path.write_text(condition)

# 5. Keep find_modules compatibility in its real implementation file.
canonical = []
for path in sorted(R_DIR.glob('*.R')):
    if path == validation_path:
        continue
    text = path.read_text()
    if re.search(r'(?m)^find_modules\.Network\s*<-\s*function\b', text):
        canonical.append(path)
if len(canonical) != 1:
    raise RuntimeError(f'expected one canonical find_modules.Network, found {canonical}')
modules_path = canonical[0]
modules = modules_path.read_text()
core_start, core_end = function_span(modules, 'find_modules.Network')
core = modules[core_start:core_end].replace(
    'find_modules.Network <- function',
    '.find_modules_network_core <- function',
    1
)
wrapper = extract_function(validation, 'find_modules.Network')
wrapper = wrapper.replace(
    '.find_modules_network_common_dictionary_base',
    '.find_modules_network_core'
)
modules = modules[:core_start] + core + '\n\n' + wrapper + modules[core_end:]
modules_path.write_text(modules)

# 6. Rename projection by responsibility and verify the fitted preprocessing identity.
projection = projection_path.read_text()
projection = projection.replace(
    '.condition_prepare_projection_input <- function(object) {',
    '.condition_prepare_projection_input <- function(object, fit) {',
    1
)
projection = projection.replace(
    'object, assay = params$rna_assay, layer = "data"',
    'object, assay = params$rna_assay, layer = fit$rna_layer',
    1
)
projection = projection.replace(
    'object, assay = params$peak_assay, layer = "data"',
    'object, assay = params$peak_assay, layer = fit$peak_layer',
    1
)
projection = replace_once(
    projection,
    '''    gene_data <- gene_data[common_cells, , drop = FALSE]
    peak_data_all <- peak_data_all[common_cells, , drop = FALSE]
    regions <- NetworkRegions(object)
''',
    '''    gene_data <- gene_data[common_cells, , drop = FALSE]
    peak_data_all <- peak_data_all[common_cells, , drop = FALSE]
    observed_fingerprint <- .condition_preprocessing_fingerprint(
        object = object, gene_data = gene_data, peak_data = peak_data_all,
        rna_layer = fit$rna_layer, peak_layer = fit$peak_layer,
        peak_value_type = fit$peak_value_type
    )
    if (!identical(observed_fingerprint, fit$preprocessing_fingerprint)) {
        stop("RNA/ATAC preprocessing identity changed after condition fitting.",
             call. = FALSE)
    }
    regions <- NetworkRegions(object)
''',
    'projection fingerprint check'
)
projection = replace_once(
    projection,
    '    prepared <- .condition_prepare_projection_input(object)\n',
    '''    required_provenance <- c(
        "rna_layer", "peak_layer", "peak_value_type",
        "preprocessing_fingerprint"
    )
    if (!all(required_provenance %in% names(fit)) ||
        any(!nzchar(vapply(fit[required_provenance], as.character,
                           character(1))))) {
        stop("The fitted condition GRN lacks preprocessing provenance.",
             call. = FALSE)
    }
    prepared <- .condition_prepare_projection_input(object, fit)
''',
    'projection fit provenance'
)
(R_DIR / 'condition_grn_projection.R').write_text(projection)
projection_path.unlink()
validation_path.unlink()

# 7. Update package description and roxygen source documentation.
description_path = ROOT / 'DESCRIPTION'
description = description_path.read_text()
description = description.replace(
    '    interaction model, the same predictor columns and the same globally\n    preprocessed RNA and ATAC layers, with scale set to FALSE.',
    '    interaction model, the same predictor columns and the same globally\n    preprocessed RNA and ATAC layers, with scale set to FALSE. The selected\n    layers, peak-value semantics and a stable preprocessing fingerprint are\n    stored in every condition fit and verified again during projection.'
)
description_path.write_text(description)

# Minimal generated-man synchronization for the newly public arguments.
for path_name in ('man/discover_grn_edges.Rd', 'man/fit_grn_from_edges.Rd',
                  'man/infer_condition_grn.Rd'):
    path = ROOT / path_name
    if not path.exists():
        continue
    text = path.read_text()
    if '\\item{rna_layer}' not in text:
        marker = '\\item{peak_to_gene_domains}'
        pos = text.find(marker)
        if pos >= 0:
            end = text.find('\n\n', pos)
            if end < 0:
                end = pos
            addition = ('\n\\item{rna_layer}{RNA assay layer used by all candidate '
                        'and fixed-dictionary fits.}\n\n'
                        '\\item{peak_layer}{ATAC assay layer used by all candidate '
                        'and fixed-dictionary fits.}\n\n'
                        '\\item{peak_value_type}{Recorded ATAC value semantics; '
                        'probability mode enforces values in [0, 1].}\n')
            text = text[:end] + addition + text[end:]
    path.write_text(text)

# 8. Add focused regression tests and permanent source guards.
test_path = ROOT / 'tests/testthat/test-common-dictionary-direct-source.R'
test_path.write_text(r'''test_that("common-dictionary functions have one direct definition", {
  files <- list.files("../../R", pattern = "[.]R$", full.names = TRUE)
  text <- vapply(files, paste, collapse = "\n", FUN.VALUE = character(1))
  names(text) <- files
  pattern <- "(?m)^([A-Za-z.][A-Za-z0-9._]*)\\s*<-\\s*function\\b"
  definitions <- unlist(lapply(text, function(value) {
    matches <- gregexpr(pattern, value, perl = TRUE)
    raw <- regmatches(value, matches)[[1L]]
    if (!length(raw) || identical(raw, character(0))) return(character())
    sub("\\s*<-\\s*function.*$", "", raw)
  }), use.names = FALSE)
  expect_false(anyDuplicated(definitions) > 0L)
  expect_false(any(basename(files) %in% c(
    "zz_common_dictionary_validation.R",
    "zz_common_dictionary_projection.R",
    "zzz_fit_schema_contract.R"
  )))
})

test_that("fixed-dictionary model controls are explicit and strict", {
  expect_match(
    paste(deparse(formals(Pando::fit_grn_from_edges)$method), collapse = ""),
    "glm", fixed = TRUE
  )
  expect_identical(formals(Pando::fit_grn_from_edges)$interaction_term, ":")
  expect_identical(formals(Pando::fit_grn_from_edges)$scale, FALSE)
  expect_true(all(c("rna_layer", "peak_layer", "peak_value_type") %in%
                  names(formals(Pando::fit_grn_from_edges))))
})
''')

workflow_path = ROOT / '.github/workflows/common-dictionary-condition-grn.yaml'
if workflow_path.exists():
    workflow = workflow_path.read_text()
    guard = '''      - name: Verify direct common-dictionary source layout
        shell: python
        run: |
          from pathlib import Path
          import re
          forbidden = {
              "zz_common_dictionary_validation.R",
              "zz_common_dictionary_projection.R",
              "zzz_fit_schema_contract.R",
          }
          files = sorted(Path("R").glob("*.R"))
          present = forbidden.intersection(path.name for path in files)
          if present:
              raise RuntimeError(f"late-loaded implementations remain: {present}")
          found = {}
          pattern = re.compile(
              r"(?m)^([A-Za-z.][A-Za-z0-9._]*)\\s*<-\\s*function\\b"
          )
          for path in files:
              for name in pattern.findall(path.read_text()):
                  found.setdefault(name, []).append(str(path))
          duplicate = {name: paths for name, paths in found.items()
                       if len(paths) > 1}
          if duplicate:
              raise RuntimeError(f"duplicate function definitions: {duplicate}")
          required = [
              ".condition_preprocessing_fingerprint",
              ".condition_validate_dictionary",
              ".condition_fit_target_matrix",
              ".condition_fit_dictionary_prepared",
              ".find_modules_network_core",
              "find_modules.Network",
              "project_condition_grn_cells",
          ]
          missing = [name for name in required if name not in found]
          if missing:
              raise RuntimeError(f"missing direct implementations: {missing}")
'''
    first_steps = '    steps:\n'
    if 'Verify direct common-dictionary source layout' not in workflow:
        workflow = workflow.replace(first_steps, first_steps + guard, 1)
    workflow_path.write_text(workflow)

# Final static audit before the workflow commits anything.
all_files = sorted(R_DIR.glob('*.R'))
found = {}
pattern = re.compile(r'(?m)^([A-Za-z.][A-Za-z0-9._]*)\s*<-\s*function\b')
for path in all_files:
    for name in pattern.findall(path.read_text()):
        found.setdefault(name, []).append(str(path))
duplicate = {name: paths for name, paths in found.items() if len(paths) > 1}
if duplicate:
    raise RuntimeError(f'duplicate functions after migration: {duplicate}')
for forbidden in (
    R_DIR / 'zz_common_dictionary_validation.R',
    R_DIR / 'zz_common_dictionary_projection.R',
    R_DIR / 'zzz_fit_schema_contract.R'
):
    if forbidden.exists():
        raise RuntimeError(f'late-loaded implementation remains: {forbidden}')
