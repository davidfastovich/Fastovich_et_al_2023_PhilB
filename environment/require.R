# Not tested if this script runs on it's own. It may require manual
# intervention to select the necessary servers for the packages.

install.packages("remotes")
library("remotes")
install_version("ncf", version = "1.3-2", repos = "http://cran.us.r-project.org")
install_version("usdm", version = "1.1-18", repos = "http://cran.us.r-project.org")
