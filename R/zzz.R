.onAttach = function(libname, pkgname) {
    if (interactive()) {
        
        ver = utils::packageVersion("maaslin3")
        
        cli::cli_inform(paste0("This is MaAsLin3 version {.strong ", ver, "}"), 
                        class = "packageStartupMessage")
        
        cli::cli_inform("{cli::symbol$bullet} {.strong Get help}: Visit the biobakery help forum at {.url https://forum.biobakery.org/}", 
                        class = "packageStartupMessage")
        
    } 
}
