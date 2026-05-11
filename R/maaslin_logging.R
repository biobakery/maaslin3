maaslin_log_levels <- c(
    NOTSET = 0,
    FINEST = 1,
    FINER = 4,
    FINE = 7,
    DEBUG = 10,
    INFO = 20,
    WARNING = 30,
    WARN = 30,
    ERROR = 40,
    CRITICAL = 50,
    FATAL = 50
)

maaslin_log_state <- new.env(parent = emptyenv())

maaslin_log_level <- function(level) {
    if (is.numeric(level)) {
        match_idx <- which(maaslin_log_levels == level)
    } else if (is.character(level)) {
        match_idx <- which(names(maaslin_log_levels) == toupper(level))
    } else {
        match_idx <- integer()
    }

    if (length(match_idx) == 0) {
        return(unname(maaslin_log_levels[["NOTSET"]]))
    }

    unname(maaslin_log_levels[match_idx[1]])
}

maaslin_log_level_name <- function(level) {
    match_idx <- which(maaslin_log_levels == level)
    if (length(match_idx) == 0) {
        return("NOTSET")
    }
    names(maaslin_log_levels)[match_idx[1]]
}

maaslin_log_compose_message <- function(msg, ...) {
    args <- list(...)

    if (!is.character(msg)) {
        msg <- paste(c(msg, args), collapse = " ")
        return(msg)
    }

    if (length(args) == 0) {
        return(msg)
    }

    args <- lapply(args, function(arg) {
        if (length(arg) != 1) {
            arg <- paste(arg, collapse = ",")
        }
        arg
    })

    do.call("sprintf", c(msg, args))
}

maaslin_log_default_format <- function(record) {
    msg <- trimws(record$msg)
    paste(record$timestamp, paste(record$levelname, record$logger, msg,
                                  sep = ":"))
}

maaslin_log_write_to_console <- function(msg, handler, ...) {
    if ("dry" %in% names(list(...))) {
        return(TRUE)
    }

    cat(paste0(msg, "\n"))
}

maaslin_log_write_to_file <- function(msg, handler, ...) {
    if ("dry" %in% names(list(...))) {
        return(exists("file", envir = handler, inherits = FALSE))
    }

    cat(paste0(msg, "\n"), file = handler$file, append = TRUE)
}

maaslin_log_get_logger <- function(name = "", ...) {
    if (!identical(name, "")) {
        stop("maaslin3 only supports the root logger")
    }

    maaslin_log_state
}

maaslin_log_get_handler <- function(handler, logger = "") {
    logger <- maaslin_log_get_logger(logger)
    if (!is.character(handler)) {
        handler <- deparse(substitute(handler))
    }

    logger$handlers[[handler]]
}

maaslin_log_remove_handler <- function(handler, logger = "") {
    logger <- maaslin_log_get_logger(logger)
    if (!is.character(handler)) {
        handler <- deparse(substitute(handler))
    }

    logger$handlers <- logger$handlers[!(names(logger$handlers) == handler)]
    invisible()
}

maaslin_log_add_handler <- function(handler, ..., logger = "") {
    logger <- maaslin_log_get_logger(logger)
    params <- list(...)

    param_names <- names(params)

    if (is.character(handler)) {
        handler_name <- handler
        if ("action" %in% names(params)) {
            action <- params[["action"]]
            params[["action"]] <- NULL
        } else if (length(params) > 0 &&
                   (is.null(param_names) || !nzchar(param_names[1]))) {
            action <- params[[1]]
            params[[1]] <- NULL
        } else {
            stop("No action for the handler provided")
        }
    } else {
        handler_name <- deparse(substitute(handler))
        action <- handler
    }

    level <- params[["level"]]
    if (is.null(level)) {
        level <- maaslin_log_levels[["NOTSET"]]
    }
    params[["level"]] <- NULL

    formatter <- params[["formatter"]]
    if (is.null(formatter)) {
        formatter <- maaslin_log_default_format
    }
    params[["formatter"]] <- NULL

    handler_env <- list2env(params, parent = emptyenv())
    handler_env$action <- action
    handler_env$level <- maaslin_log_level(level)
    handler_env$formatter <- formatter

    maaslin_log_remove_handler(handler_name)
    if (isTRUE(action(NA, handler_env, dry = TRUE))) {
        logger$handlers[[handler_name]] <- handler_env
    }

    invisible(handler_env)
}

maaslin_log_set_level <- function(level, container = "") {
    if (is.null(container)) {
        stop("NULL container provided: cannot set level for NULL container")
    }

    if (is.character(container)) {
        container <- maaslin_log_get_logger(container)
    }

    container$level <- maaslin_log_level(level)
    invisible()
}

maaslin_log_basic_config <- function(level = 20) {
    maaslin_log_set_level(level)
    maaslin_log_add_handler("basic.stdout", maaslin_log_write_to_console)
    invisible()
}

maaslin_log_reset <- function() {
    maaslin_log_state$name <- ""
    maaslin_log_state$level <- maaslin_log_level("INFO")
    maaslin_log_state$handlers <- list()
    invisible()
}

maaslin_log_record <- function(level, msg, ..., logger = "") {
    logger_env <- maaslin_log_get_logger(logger)
    numeric_level <- maaslin_log_level(level)

    if (numeric_level < logger_env$level) {
        return(invisible(FALSE))
    }

    record <- list(
        msg = maaslin_log_compose_message(msg, ...),
        timestamp = Sys.time(),
        logger = logger,
        level = numeric_level,
        levelname = maaslin_log_level_name(numeric_level)
    )

    for (handler in logger_env$handlers) {
        if (numeric_level >= handler$level) {
            handler$action(handler$formatter(record), handler, record)
        }
    }

    invisible(TRUE)
}

maaslin_logdebug <- function(msg, ..., logger = "") {
    maaslin_log_record("DEBUG", msg, ..., logger = logger)
}

maaslin_loginfo <- function(msg, ..., logger = "") {
    maaslin_log_record("INFO", msg, ..., logger = logger)
}

maaslin_logwarn <- function(msg, ..., logger = "") {
    maaslin_log_record("WARN", msg, ..., logger = logger)
}

maaslin_logerror <- function(msg, ..., logger = "") {
    maaslin_log_record("ERROR", msg, ..., logger = logger)
}

maaslin_log_reset()
maaslin_log_basic_config()
