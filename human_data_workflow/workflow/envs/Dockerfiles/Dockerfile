from r-base:3.6.1

run apt-get update && \
  apt-get install -y libcurl4-openssl-dev libssl-dev libssh2-1-dev libxml2-dev zlib1g-dev

run R -e "install.packages(c('devtools'))"

run R -e "install.packages('BiocManager')"
run R -e "BiocManager::install('multtest')"

run R -e "devtools::install_url('https://cran.r-project.org/src/contrib/Archive/XML/XML_3.98-1.19.tar.gz')"

run R -e "BiocManager::install('annotate')"
run R -e "BiocManager::install('GSEABase')"
run R -e "BiocManager::install('GSVA')"

run R -e "devtools::install_url('https://cran.r-project.org/src/contrib/Archive/mnormt/mnormt_1.5-6.tar.gz')"
run R -e "devtools::install_url('https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_3.1.3.tar.gz')"
run R -e "devtools::install_github('dviraran/SingleR', ref = '195ffc24b213c6b3b11395f86e88ed60238dfbc4')"

