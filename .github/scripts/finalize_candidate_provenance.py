from pathlib import Path
import re

path = Path('R/condition_grn.R')
text = path.read_text()


def function_span(source, name):
    match = re.search(r'(?m)^' + re.escape(name) + r'\s*<-\s*function\b', source)
    if match is None:
        raise RuntimeError(f'missing function: {name}')
    if re.search(r'(?m)^' + re.escape(name) + r'\s*<-\s*function\b', source[match.end():]):
        raise RuntimeError(f'duplicate function before patch: {name}')
    opening = source.find('{', match.end())
    depth = 0
    quote = None
    escaped = False
    comment = False
    for index in range(opening, len(source)):
        char = source[index]
        if comment:
            if char == '\n':
                comment = False
            continue
        if quote is not None:
            if escaped:
                escaped = False
            elif char == '\\':
                escaped = True
            elif char == quote:
                quote = None
            continue
        if char == '#':
            comment = True
        elif char in ('"', "'"):
            quote = char
        elif char == '{':
            depth += 1
        elif char == '}':
            depth -= 1
            if depth == 0:
                return match.start(), index + 1
    raise RuntimeError(f'unclosed function: {name}')


def replace_once(source, old, new, label):
    count = source.count(old)
    if count != 1:
        raise RuntimeError(f'{label}: expected one match, found {count}')
    return source.replace(old, new, 1)

# Candidate discovery records the exact shared numerical representation.
start, end = function_span(text, '.condition_discover_edges_prepared')
fun = text[start:end]
old = '''    attr(out, "source_label") <- source_label
    attr(out, "source_type") <- source_type
    out
'''
new = '''    attr(out, "source_label") <- source_label
    attr(out, "source_type") <- source_type
    attr(out, "rna_layer") <- prepared$rna_layer
    attr(out, "peak_layer") <- prepared$peak_layer
    attr(out, "peak_value_type") <- prepared$peak_value_type
    attr(out, "preprocessing_fingerprint") <-
        prepared$preprocessing_fingerprint
    attr(out, "dictionary_input_schema") <-
        "pando_candidate_input_provenance_v1"
    out
'''
fun = replace_once(fun, old, new, 'candidate provenance')
text = text[:start] + fun + text[end:]

# Exact union rejects mixed preprocessing references and carries provenance.
start, end = function_span(text, 'union_grn_edges')
fun = text[start:end]
old = '''    all_rows <- do.call(rbind, lapply(seq_along(inputs), function(i) {
'''
new = '''    provenance_fields <- c(
        "rna_layer", "peak_layer", "peak_value_type",
        "preprocessing_fingerprint", "dictionary_input_schema"
    )
    provenance <- lapply(inputs, function(candidate) {
        stats::setNames(lapply(provenance_fields, function(field) {
            attr(candidate, field, exact = TRUE)
        }), provenance_fields)
    })
    complete_provenance <- vapply(provenance, function(value) {
        all(vapply(value, function(item) {
            is.character(item) && length(item) == 1L &&
                !is.na(item) && nzchar(item)
        }, logical(1)))
    }, logical(1))
    if (any(complete_provenance) && !all(complete_provenance)) {
        stop(
            "Candidate tables mix verified and unverified preprocessing ",
            "provenance.", call. = FALSE
        )
    }
    if (all(complete_provenance)) {
        reference <- provenance[[1L]]
        same_reference <- vapply(provenance, function(value) {
            identical(value, reference)
        }, logical(1))
        if (!all(same_reference)) {
            stop(
                "Global and condition candidates use different RNA/ATAC ",
                "layers, value semantics, or preprocessing fingerprints.",
                call. = FALSE
            )
        }
    }

    all_rows <- do.call(rbind, lapply(seq_along(inputs), function(i) {
'''
fun = replace_once(fun, old, new, 'union provenance validation')
old = '''    class(dictionary) <- c("PandoEdgeDictionary", "data.frame")
    dictionary
'''
new = '''    class(dictionary) <- c("PandoEdgeDictionary", "data.frame")
    attr(dictionary, "preprocessing_provenance_verified") <-
        all(complete_provenance)
    if (all(complete_provenance)) {
        for (field in provenance_fields) {
            attr(dictionary, field) <- provenance[[1L]][[field]]
        }
    }
    dictionary
'''
fun = replace_once(fun, old, new, 'union provenance output')
text = text[:start] + fun + text[end:]

# Fixed fitting rejects a dictionary made in another numerical reference.
start, end = function_span(text, '.condition_validate_dictionary')
fun = text[start:end]
old = '''    expected <- paste(
        dictionary$target, dictionary$tf, dictionary$region, sep = "||"
    )
'''
new = '''    dictionary_fingerprint <- attr(
        dictionary, "preprocessing_fingerprint", exact = TRUE
    )
    dictionary_rna_layer <- attr(dictionary, "rna_layer", exact = TRUE)
    dictionary_peak_layer <- attr(dictionary, "peak_layer", exact = TRUE)
    dictionary_peak_value_type <- attr(
        dictionary, "peak_value_type", exact = TRUE
    )
    provenance_values <- list(
        dictionary_fingerprint, dictionary_rna_layer,
        dictionary_peak_layer, dictionary_peak_value_type
    )
    has_provenance <- vapply(provenance_values, function(value) {
        is.character(value) && length(value) == 1L &&
            !is.na(value) && nzchar(value)
    }, logical(1))
    if (any(has_provenance) && !all(has_provenance)) {
        stop("The dictionary contains incomplete preprocessing provenance.",
             call. = FALSE)
    }
    if (all(has_provenance) &&
        (!identical(dictionary_fingerprint,
                    prepared$preprocessing_fingerprint) ||
         !identical(dictionary_rna_layer, prepared$rna_layer) ||
         !identical(dictionary_peak_layer, prepared$peak_layer) ||
         !identical(dictionary_peak_value_type,
                    prepared$peak_value_type))) {
        stop(
            "The frozen dictionary and fixed fit use different RNA/ATAC ",
            "preprocessing references.", call. = FALSE
        )
    }

    expected <- paste(
        dictionary$target, dictionary$tf, dictionary$region, sep = "||"
    )
'''
fun = replace_once(fun, old, new, 'fit provenance validation')
text = text[:start] + fun + text[end:]

# Persist dictionary-level verification in every standard Network.
old = '''            preprocessing_fingerprint = prepared$preprocessing_fingerprint,
            adjust_method = adjust_method,
'''
new = '''            preprocessing_fingerprint = prepared$preprocessing_fingerprint,
            dictionary_preprocessing_provenance_verified = isTRUE(attr(
                edge_dictionary,
                "preprocessing_provenance_verified",
                exact = TRUE
            )),
            adjust_method = adjust_method,
'''
text = replace_once(text, old, new, 'Network dictionary provenance')

path.write_text(text)

# Extend focused regression coverage without requiring a GRNData fixture.
test_path = Path('tests/testthat/test-common-dictionary-direct-source.R')
test = test_path.read_text()
addition = r'''

test_that("exact edge union enforces one preprocessing reference", {
  candidate <- data.frame(
    target = "G", tf = "TF", region = "chr1-1-2",
    atac_feature_id = "chr1-1-2", peak_target_cor = 0.2,
    tf_target_cor = 0.3, source_label = "x",
    source_type = "condition", edge_id = "G||TF||chr1-1-2",
    stringsAsFactors = FALSE
  )
  class(candidate) <- c("PandoEdgeDictionary", "data.frame")
  stamp <- function(value, fingerprint) {
    attr(value, "rna_layer") <- "data"
    attr(value, "peak_layer") <- "data"
    attr(value, "peak_value_type") <- "normalized"
    attr(value, "preprocessing_fingerprint") <- fingerprint
    attr(value, "dictionary_input_schema") <-
      "pando_candidate_input_provenance_v1"
    value
  }
  global <- stamp(candidate, "same")
  condition <- stamp(candidate, "same")
  union <- Pando::union_grn_edges(
    global, list(control = condition)
  )
  expect_true(isTRUE(attr(
    union, "preprocessing_provenance_verified", exact = TRUE
  )))
  expect_identical(
    attr(union, "preprocessing_fingerprint", exact = TRUE), "same"
  )

  mismatched <- stamp(candidate, "different")
  expect_error(
    Pando::union_grn_edges(global, list(control = mismatched)),
    "different RNA/ATAC"
  )
  unverified <- candidate
  attributes(unverified)[c(
    "rna_layer", "peak_layer", "peak_value_type",
    "preprocessing_fingerprint", "dictionary_input_schema"
  )] <- NULL
  expect_error(
    Pando::union_grn_edges(global, list(control = unverified)),
    "mix verified and unverified"
  )
})
'''
if 'exact edge union enforces one preprocessing reference' not in test:
    test_path.write_text(test.rstrip() + addition + '\n')

# Static final audit.
source = Path('R/condition_grn.R').read_text()
for required in (
    'pando_candidate_input_provenance_v1',
    'preprocessing_provenance_verified',
    'different RNA/ATAC ',
    'layers, value semantics, or preprocessing fingerprints.',
    'frozen dictionary and fixed fit use different RNA/ATAC '
):
    if required not in source:
        raise RuntimeError(f'missing candidate provenance contract: {required}')
