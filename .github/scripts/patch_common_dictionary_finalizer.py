from pathlib import Path

path = Path('.github/scripts/finalize_common_dictionary_contract.py')
text = path.read_text()
old = '''condition = replace_once(
    condition,
    '    peak_to_gene_method <- match.arg(peak_to_gene_method)\\n',
    ''' + "'''" + '''    peak_to_gene_method <- match.arg(peak_to_gene_method)
    peak_value_type <- match.arg(peak_value_type)
    valid_layer <- function(value) {
        is.character(value) && length(value) == 1L && !is.na(value) &&
            nzchar(trimws(value))
    }
    if (!valid_layer(rna_layer) || !valid_layer(peak_layer)) {
        stop("RNA and ATAC layer names must be non-empty strings.",
             call. = FALSE)
    }
''' + "'''" + ''',
    'layer validation'
)
prepare_start, prepare_end = function_span(condition, '.condition_prepare_common_input')
prepare = condition[prepare_start:prepare_end]
'''
new = '''prepare_start, prepare_end = function_span(condition, '.condition_prepare_common_input')
prepare = condition[prepare_start:prepare_end]
prepare = replace_once(
    prepare,
    '    peak_to_gene_method <- match.arg(peak_to_gene_method)\\n',
    ''' + "'''" + '''    peak_to_gene_method <- match.arg(peak_to_gene_method)
    peak_value_type <- match.arg(peak_value_type)
    valid_layer <- function(value) {
        is.character(value) && length(value) == 1L && !is.na(value) &&
            nzchar(trimws(value))
    }
    if (!valid_layer(rna_layer) || !valid_layer(peak_layer)) {
        stop("RNA and ATAC layer names must be non-empty strings.",
             call. = FALSE)
    }
''' + "'''" + ''',
    'layer validation'
)
'''
if text.count(old) != 1:
    raise RuntimeError(f'expected one finalizer block, found {text.count(old)}')
path.write_text(text.replace(old, new, 1))
