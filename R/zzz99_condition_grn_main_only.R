# Remove stability and sensitivity refits from the canonical runtime path.

.condition_fit_targets_reference <- .condition_fit_targets

.condition_refit_stability_diagnostics <- function(...) NULL

.condition_fit_targets <- function(...) {
    answer <- .condition_fit_targets_reference(...)
    answer$fits <- lapply(answer$fits, function(result) {
        contract <- result$fit_contract
        if (is.list(contract)) {
            contract$refit_stability <- NULL
            if (is.list(contract$refit)) {
                contract$refit$stability_status <- NULL
            }
            result$fit_contract <- contract
        }
        result
    })
    answer
}
