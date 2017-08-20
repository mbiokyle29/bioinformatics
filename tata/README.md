# TATA - Finding TATA* boxes with CAGEseq data

## Pre-reqs

A BSgenome object of your ref seq. If you don't have one already:

1. Create a seed file (see BSgenome)
2. Fire up R

    ```R
        library(BSGenome)
        forgeBSgenomeDataPkg("path/to/my/seed")
    ```
3.  Then build the package (shell)
    
    ```bash
        R CMD build <package made above>
        R CMD check <tar of package>
        R CMD INSTALL <tar of package> <dest of R libs>

    ```
4. I had a bit of trouble with that approach. Opening up R in doing install.packages("./<tarball>")

## Analysis

1. CAGEr - Use CAGEr identify the actual TSS's from bam file