# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE


.onAttach <- function(libname, pkgname) {
  if (require("R.oo")) {
    pkg <- Package(pkgname);
    assign(pkgname, pkg, pos=getPosition(pkg));
  }

  pi <- utils::packageDescription(pkgname);
  packageStartupMessage(pkgname, " v", pi$Version, " (", 
    pi$Date, ") successfully loaded. See ?", pkgname, " for help."); 
}
