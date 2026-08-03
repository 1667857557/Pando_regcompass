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
text = text.replace(old, new, 1)

old_projection_write = '''(R_DIR / 'condition_grn_projection.R').write_text(projection)
projection_path.unlink()
validation_path.unlink()
'''
new_projection_write = '''# Projection is implemented only in its responsibility-named source file.
project_start, project_end = function_span(
    condition, 'project_condition_grn_cells'
)
condition = condition[:project_start] + condition[project_end:]
condition_path.write_text(condition)
(R_DIR / 'condition_grn_projection.R').write_text(projection)
projection_path.unlink()
validation_path.unlink()
'''
if text.count(old_projection_write) != 1:
    raise RuntimeError(
        f'expected one projection write block, found '
        f'{text.count(old_projection_write)}'
    )
text = text.replace(old_projection_write, new_projection_write, 1)
path.write_text(text)
