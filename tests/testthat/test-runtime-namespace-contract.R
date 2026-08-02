namespace_references <- function(files) {
    references <- list()
    add_reference <- function(operator, package, symbol, file) {
        references[[length(references) + 1L]] <<- data.frame(
            operator = operator,
            package = package,
            symbol = symbol,
            file = basename(file),
            stringsAsFactors = FALSE
        )
    }
    walk <- function(value, file) {
        if (is.call(value)) {
            head <- value[[1L]]
            operator <- if (is.symbol(head)) as.character(head) else ""
            if (operator %in% c("::", ":::") && length(value) >= 3L) {
                add_reference(
                    operator,
                    as.character(value[[2L]]),
                    as.character(value[[3L]]),
                    file
                )
            }
            invisible(lapply(as.list(value), walk, file = file))
        } else if (is.expression(value) || is.pairlist(value)) {
            invisible(lapply(as.list(value), walk, file = file))
        }
        invisible(NULL)
    }
    for (file in files) {
        walk(parse(file = file, keep.source = FALSE), file)
    }
    if (!length(references)) {
        return(data.frame(
            operator = character(), package = character(),
            symbol = character(), file = character()
        ))
    }
    unique(do.call(rbind, references))
}

split_description_dependencies <- function(value) {
    if (!length(value) || is.na(value) || !nzchar(trimws(value))) {
        return(character())
    }
    entries <- trimws(unlist(strsplit(value, ",", fixed = TRUE)))
    entries <- trimws(sub("\\s*\\(.*\\)$", "", entries))
    entries[nzchar(entries)]
}

namespace_import_packages <- function(lines) {
    pattern <- "^\\s*import(?:From|ClassesFrom)?\\(([^,)]+)"
    matches <- regexec(pattern, lines, perl = TRUE)
    values <- regmatches(lines, matches)
    unique(vapply(
        values[lengths(values) >= 2L],
        function(value) gsub("['\"]", "", trimws(value[[2L]])),
        character(1)
    ))
}

test_that("runtime namespace dependencies are direct package dependencies", {
    root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    description <- read.dcf(file.path(root, "DESCRIPTION"))[1L, , drop = TRUE]
    declared <- unique(c(
        split_description_dependencies(description[["Depends"]]),
        split_description_dependencies(description[["Imports"]])
    ))
    declared <- setdiff(declared, "R")

    namespace_packages <- namespace_import_packages(
        readLines(file.path(root, "NAMESPACE"), warn = FALSE)
    )
    source_files <- list.files(
        file.path(root, "R"), pattern = "\\.[Rr]$",
        recursive = TRUE, full.names = TRUE
    )
    references <- namespace_references(source_files)
    runtime_packages <- unique(c(namespace_packages, references$package))
    missing <- setdiff(runtime_packages, c(declared, "base", "Pando"))

    expect_length(
        missing,
        0L,
        info = paste(
            "Runtime namespaces missing from Depends/Imports:",
            paste(sort(missing), collapse = ", ")
        )
    )
})

test_that("runtime double-colon calls address exported objects", {
    root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    source_files <- list.files(
        file.path(root, "R"), pattern = "\\.[Rr]$",
        recursive = TRUE, full.names = TRUE
    )
    references <- namespace_references(source_files)
    public <- references[references$operator == "::", , drop = FALSE]

    failures <- character()
    for (index in seq_len(nrow(public))) {
        package <- public$package[[index]]
        symbol <- public$symbol[[index]]
        if (!requireNamespace(package, quietly = TRUE)) {
            failures <- c(
                failures,
                paste0(package, "::", symbol, " (namespace unavailable)")
            )
        } else if (!symbol %in% getNamespaceExports(package)) {
            failures <- c(
                failures,
                paste0(
                    package, "::", symbol, " in ", public$file[[index]],
                    " (object is not exported)"
                )
            )
        }
    }

    expect_length(
        unique(failures),
        0L,
        info = paste(unique(failures), collapse = "; ")
    )
})

test_that("runtime code does not reach into dependency internals", {
    root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    source_files <- list.files(
        file.path(root, "R"), pattern = "\\.[Rr]$",
        recursive = TRUE, full.names = TRUE
    )
    references <- namespace_references(source_files)
    internal <- references[references$operator == ":::", , drop = FALSE]

    expect_equal(
        nrow(internal),
        0L,
        info = paste(
            paste0(
                internal$package, "::: ", internal$symbol,
                " in ", internal$file
            ),
            collapse = "; "
        )
    )
})
