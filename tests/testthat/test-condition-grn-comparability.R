test_that('reference comparison mask separates zero from not estimable', {
    eligibility <- matrix(
        c(
            TRUE, TRUE,
            TRUE, FALSE,
            FALSE, TRUE
        ),
        nrow = 3L,
        byrow = TRUE,
        dimnames = list(
            c('edge_both', 'edge_reference_only', 'edge_condition_only'),
            c('Control', 'Drug')
        )
    )

    observed <- Pando:::.condition_reference_comparison_mask(
        eligibility, 'Control'
    )

    expect_true(observed['edge_both', 'Drug'])
    expect_false(observed['edge_reference_only', 'Drug'])
    expect_false(observed['edge_condition_only', 'Drug'])
    expect_true(all(observed[, 'Control'] == eligibility[, 'Control']))
    expect_identical(dimnames(observed), dimnames(eligibility))
})

test_that('comparison mask validates the reference condition', {
    eligibility <- matrix(
        TRUE, nrow = 1L, ncol = 2L,
        dimnames = list('edge1', c('Control', 'Drug'))
    )

    expect_error(
        Pando:::.condition_reference_comparison_mask(eligibility, 'Missing'),
        'reference_condition'
    )
})

test_that('combined fit contract records comparison semantics', {
    body_text <- paste(
        deparse(body(Pando:::.condition_combine_fit_contracts)),
        collapse = '\n'
    )
    expect_match(body_text, 'comparison_mask', fixed = TRUE)
    expect_match(body_text, 'eligibility_mask', fixed = TRUE)
})
