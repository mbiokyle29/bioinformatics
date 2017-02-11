# Bioinformatics Manual

- This manual represents a summary of how I would approach doing bioinformatics programming and analysis. It is basically a collection of the things that I have found which work best in my several years as a DIY bioinformatics analyst.

- A huge **DISCLAIMER** first. The biggest thing to stress is that nothing in here is a requirement nor it is garunteed to always work. At this point general bioinformatics best practices don't really exist, so take everything with a grain of salt. There is no perfect way to run the Tophat aligner, no perfect set of parameters to get the right alignment from BWA, _etc_, _etc_.

- I absolutely learned a lot of this from other people, hopefully they are not mad that I don't name tyh

### Chapter #1 - General Strategies / Good Ideas

###### I am first going to highlight a few general things to keep in mind while working on analysis projects. The thing to always keep in your mind is this:
#
#

>If I were to stop working on this project for 3 weeks, and come back to it would I know what was in the projects directory?

I cannot count how many times this has happened to me. You return after working on something else and don't remember what anything is. Here is my approach.
1. Always have an `ANALYSIS.md` file at the root of the project directory. Keep good notes of what you have done, what the project is about and anything else important. One big thing here is tracking any renaming events. When you get new sequence data, the filenames will often be not helpful. Record a mapping of the original filenames to the new names (and include the actual commands you ran to change the names)
2. Develop a consistent directory structure and use it for every project!
    -  Here is a break down of what I use:
       - `fastq/` - Where I keep the raw sequence data
       - `alignments/` - Where I store results from the alignment stage a of a project, usually with RNAseq I have bams and RSEM count results, alignment statistics, etc
       - `docs/` - Where I put anything document like, so things like pdf, and excel tables. I loathe excel but people don't love getting plain text files so that's that.
       - `R/` - I usually have this to store all the R code I write. Usually I have multiple subdirectories under `R/` that cover various analysis approachs or what have you.
       - `scripts/` - This is for bash and other non-R scripts for any kind of one-off things or what have you. I like to save any renaming commands into a script and save it here just so it's around.
       - `var/` - This is the trashbin where I put old stuff, or no longer needed stuff. I try not to delete anything, it goes in here.
3. The protection of important data is usually your primary concern. I have on multiple occasionas ran an `rm -f` and nuked a bunch of fastq files, which is a costly mistake. Due to our somewhat limited resources, we cannot nessecarily back up everything. So a few tips:
    - Establish a set of absolutley nessecary files and back those up. In 99% of the cases these will be just the fastq files. As long as you document your analysis work, you can and should be able to recreate things from just fastq files. So if you delete `.bam` files by accident you can just realign.
    - Use linux file permissions to protect important things. You can/should make `fastq` files readonly, so that you dont accidently modify them. Realistically there should be no need to modify `fastq` files: `chmod a-w *.fastq` should make the file unwritable (assuming my chmod syntax is correct).
    - If you want to reall go crazy with it you can do `sudo chattr +i /path/to/my/file.fastq`. This sets the immutable flag on the file and keeps it from being changed in anyway even by root.
    - Try to never move fastq files around! Alot of the tools I have written and some other ones work with the idea of: *Give me a directory and I'll process all the `.fastq` files in it*. This works well in most cases but sometimes we want to mess around or just align a subset of files. So we can create a new directory, symlink the fastqs we want a go from there.
        ```bash
        mkdir fastq-subset/
        cd fastq-subset/
        
        # repeat this step for all the fastqs you want to use
        # note that you should always fullpaths
        ln -s /full/path/to/fastq/file.fastq /full/path/to/fastq-subset/file.fastq
        ```
4. Always try and visualize things as much as possible. Whether that means looking directly at the alignment step results using a bam/wig viewer or generating some plots in R or excel, too it as often as possible. Many **A Ha** moments have happend from just stoping, taking a step back and lookin at intermediate steps in the analysis. The exact details about how to do this depends on the situation. but simple XY plots are often useful.

5. Naming conventions are a nice tool to keep things organized and allow files to be easily identififed. For any kind of "data" file we break it down as such:
    - `<Assay-Type>_<Cell Line / type>_<Treatment>-rep#.file_type`
    - Assay type specifies whether this is `DNAseq`, `RNAseq`, `CHiPSeq` etc.
    - Cell Line / Type describes the cells, I usually have pick a useful term (assuming you understand the project) as opposed to the actual cell line but sometimes the actual cell.
        - Examples: `noks`, `noks-akata`, `p14p16KO`
    - Treatment specifies any drug treatment or other modification made to the cell, use this for drugs and targetted mutations, etc. If you have a drug treatment with plus/minus specify here (i.e. don't have leave out the drug treatment in the name of the negative condition). If you have two treatments or more and a negative condition define it as `no_treatment` or something.
        - Examples: `CaFBS`, `mIRF4`, `4HT-pos`, `NoTrt`
    - Finally use standard file type extensions. Don't add a `.txt` at the end of these and don't use crappy text editors that do that either. I prefer:
        - fastq --> `.fastq`
        - fasta --> `.fasta`
        - multi-fasta --> `.mfasta`
        - binary alignment --> `.bam`
        - alignment --> `.sam`
        - bed --> `.bed`
        - wig --> `.wig`
        - bigWig --> `.bw`
        - bedGraph --> `.bg`
6. Software Versioning! Often times I will run into an issue with a bioformatics tool, after being a good citizen and debugging for a few hours I will often seek out assistance on various help sites, either general boards like SO/Biostars or directly to the projects help page (eg: Google Group). The #1 response I get is what version are you using? I am a lazy person and will often do this:
    ```shell
    sudo apt-get install bowtie
    ```
    over this:
    ```
    wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip
    unzip bowtie-1.1.2-linux-x86_64.zip
    cd bowtie-1.1.2
    # blah
    # blah
    ```
    But the second option is almost always the best, since `apt` packages are usually woefully out of date. Here is my strategy to force my self to use option 2:
    - Pick a central location to install applications, one that is not in your home dir so it can be accessiable to other users (this depends largely on your setup I suppose)
    - Install packages into the dir, and maintain their version numbers. So let's say I install `bowtie-1.1.2` into my application directory (for me: `/data/app/`).
        ```bash
        ls /data/app/
        bowtie-1.1.2/
        ```
        Since we just downloaded this, we know it's the latest version so we can make a symlink for it
        ```bash
        ln -s /data/app/bowtie-1.1.2/ /data/app/bowtie
        ```
        Now `/data/app/bowtie/` will point to our most up to date version of bowtie and we don't have to worry about which version it actually is. When we get a new version we just update the symlink. Then we can also symlink the executable `/data/app/bowtie/bowtie` in `/usr/bin/` or somewhere else convenient and it all should just work.
7. Using the right tool for the job! This one is huge. The following examples should suffice. At one point we were working on rebuilding an EBV genome. Our actual EBV strain had a GFP (florecent) gene inserted into the middle. So we were using BLAST and a few other tools to "map in" the GFP to our reference. We were passing it back and forth. It was being transferred back and forth as a `.word` document. Since I was on linux I was opening it in LibreOffice, copying everything into Sublime Text (plain text!) and going from there. Well LibreOffice has a **silent** line length limit which meant we were losing pieces of the reference as we passed it around and it took quite some time to realize it was happening and to track it down. 

    There are a number of human genes with names like SEPT1, OCT10, etc. When you paste those into Excel it will sometimes format them into dates like `9/1` and `10/10`. This is actually a very widespread problem: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-80
    - So: USE A PLAIN TEXT EDITOR ALL THE TIME NO MATTER WHAT
    - You will eventually need to dump stuff to excel, but please try to minimize that. Find a feature rich plain text editor and use it. When your slicing `.tsv` files and writing python code you need the range of support you get from plain text and I would not reccomend an IDE. I will make the case for sublime text, it has great find and replace with regex multiline and column selection and a number of installable packages for bioinformatics (the nucleotide syntax highlighting is the best). If you don't want Sublime and are on `OSX` TextWrangler is pretty nice, if you are for some reason on windows Notepad++ is your best bet.

     - Use a good Terminal / Shell / Terminal Emulator
    And master it! I use `rxvt-unicode` which is very light weight, I can spawn new tabs with `Shift-Down` and move between tabs with `Shift-<Left/Right>`. There are a few other bells and whistles but not much. I rely heavily on custom bash aliases and functions. I would figure out how to make these and build up your own set.

        Here are a few good ones (most of these were given to me by various ppl):
        ```bash
        # nice catting of .tsv files
        function tsv() { cat "$@" | column -t -s"\t" | less -S; }
        
        # Colorize ls
        alias ls="ls --color=auto"
        
        # "Super" Search
        alias lss="ls -tlhga --color=auto"
        
        # clear the screen
        alias c="clear"
    
        # up a dir
        alias ..="cd .."
### Chapter #2 - Doing an Analysis Project
###### Where I outline what a typical project looks like, touch on a few basics of bioinformatics files and tools you will likely use along the way.
#
The first step of any project will usually be `.fastq` files. A quick reminded, `.fastq` files contain raw sequence data from a sequence machine. It will 99% of the time be short reads (50-100) base pair reads. Each read is represented as a four line record with the following format:
1. A header with info about the seqeuence machine, lane, run and barcode.
2. The acutal sequenced bases
3. The most wasteful line in all of bioinformatics
4. The Q scores for each base calls, **usually** using the following scoring (left is best, right is worst)
> !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
#

```
@MACHINE_NAME ###  ###:### AATGCATCAG
ACTAGCTAGCTAGCTGATCGATCGATGCTAGCTAGCTAGCTAG
+ (sometimes info from line 1 is here again)
Q scores, same length as read
```

If you get the `.fastq` files from the UW biotech center, they will include `FastQC` plots for each sample. You should look at these results to make sure the raw data looks good.

The following has good explainations of what each graph means / how to interpret:
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/

Before we run an alignment (usually the first step) we need to get our reference genome and alignment indicies built. The starting point for alignment indicies is usually a `.fasta` file. A `.fasta` file contains the sequence of an entire genome (sometimes broken down by chromosome). A fasta file looks like this:
```
> Header with information about the sequence
ATCATGCTAGCTAGCTTTACGATCG.......
````

Basically most alignment programs require you to prebuild data files for an alignment. To make a bowtie 1/2 index run the following command:
```bash
bowtie(2)-build $reference_fasta $basename
```
Note that basename can be a path ending in a "basename". If we passed `/data/refs/organism-name/genome-name` as the basename, it would make all the index files in `/data/refs/organism-name/` and each file would be prepended with `genome-name`. I always do this in `/data/refs/` and make a subdirectory for each organism/build. Let run a concrete example.
```bash
# pretend we have a new human genome hg20 we want to work with
cd /data/refs/
mkdir hg20
cd hg20
mv ~/Downloads/hg20.fasta .
bowtie-build hg20.fasta hg20
```
Now all the index files live in `/data/refs/hg20/` and we are ready to run the aligner. I am showing these examples as all bare minimum options, you will need to consult the manuals and the project specifications to determine the correct (best guess) set of flags.
```bash
cd /data/projects/hg20-project/
# the format is: bowtie $index $fastq-file $outfile
# note we give the path+basename as the index
bowtie /data/refs/hg20/hg20 /data/projects/hg20-project/fastq/DNAseq-noks-rep1.fastq /data/projects/alignments/DNAseq-noks-rep1.sam
```

Let's say we have a while bunch of fastq files that we want to align, we could make a bash script with the alignment command written out for each file but that is a pain. Since we have good naming conventions we can use bash wildcards to make it easier. We will use a loop and also a string substitution to make the name of the sam file.
```bash
for file in /data/projects/hg20-project/fastq/*.fastq
    do
        bowtie /data/refs/hg20/hg20 /data/projects/hg20-project/fastq/$file /data/projects/hg20-project/alignments/${file/\.fastq/\.sam}
done;
```
Okay, so a few quick bash notes first the name of the matching file gets saved in `$file` during each iteration and bash lets us do string substitution in line using the following syntax: `${VARIABLE_NAME/MATCH_STRING/REPLACE_STRING}`. Note I escaped the dots in the file names with `\` which is probably not 100% nessecary. If you are planning on actually running something like this I would reccommend first running it using `echo` to make sure everything works right.
```bash
for file in *.fastq
    do
        # test your fastqs are being found
        echo $file
        
        # test your string rep is working
        echo ${file/\.fastq/\.sam}
done;
```
At this point we will have alignment files in either `sam` or `bam` format. `bam` is just compressed `sam` which we will want to have most of the time. `samtools` is our go to for working with `sam`/`bam` files. Assuming we have `sam` files, we usually want to convert to `bam`, sort and index.

```bash
samtools view -bS -o outfile.bam $sam_file

# samtools sort is a bit weird with output
# it will add .bam to whatever you give it
# I usually do this, to have it called .sorted.bam
# samtools adds the .bam after the .sorted
samtools sort $bam_file ${bam_file/\.bam/\.sorted}

samtools index $sorted_bam_file
```

At this point we have our alignments completed, and bam tools indexed for downstream tools. At this point it really depends on the project, but briefly:
- If this is DNAseq and we are doing variant calling we will likely run the `bam` files through GATK, 

# Chapter #3 - An R primer
###### A quick lesson in R basics and bioinformatics stuff.
#
`R` or `R lang` is a programming language and enviornment used very much in the staistics world.
You will hopefully not be writing much actual `R` code in the sense of modules/packages but mostly analysis scripts which run other packages and make some plots. The best way to do `R` work is to use an interactive `R` REPL to try stuff and work through the issues you run into and then finalize the code you use. Sublime has a great REPL, otherwise I'd use R studio. If you must, the `R` terminal console is not awful...

Most R sessions start with loading in some data from a plain text file. Most of the time the data will roughly be in a "table" like format. Think of an excel table with rows and columns, each column may or may not have a Header or title, and rows might also have labels. We will use the `R` command `read.table` function, and be sure to set the params correctly.

```R
# in the R console we can use ? to get info about a function
?read.table

read.table(file, header = FALSE, sep = "", quote = "\"'",
                dec = ".", numerals = c("allow.loss", "warn.loss", "no.loss"),
                row.names, col.names, as.is = !stringsAsFactors,
                na.strings = "NA", colClasses = NA, nrows = -1,
                skip = 0, check.names = TRUE, fill = !blank.lines.skip,
                strip.white = FALSE, blank.lines.skip = TRUE,
                comment.char = "#",
                allowEscapes = FALSE, flush = FALSE,
                stringsAsFactors = default.stringsAsFactors(),
                fileEncoding = "", encoding = "unknown", text, skipNul = FALSE)

```
We can see the `read.table` command takes 1 requried/positional argument and a whole bunch of named optional args. There is alot here but all you will really need to focus on are:
  - `file`: This is a path to the file we want to load. R does have a working directory and changing but best to just use fullpaths here.
  - `header`: This is a boolean arg which specifies if the first row is column names. `R`'s booleans are `TRUE` and `FALSE` but you can shorten them to `T` and `F` for lazyness. Set this to `T` most of the time
  - `row.names`: This allows us to pass a vector (more or less a list) of row names OR set it to a column index. In most cases our first column will be the row names, so we would do `row.names=1` (R is not a big fan of 0 index)
    - A warning about rownames, anytime you specify rownames they must be unique!
    - If you don't specify, they default to an incrementing integer count
  - `sep`: This is the field seperator which will usually be "\t" or ","
The return value should be saved to a variable using either `<-` or `=`. There are some differences that I don't remember but they should be more or less interchangable. The `<-` is more fun though!

```R
# pure data file no row names or headers
# R ppl like to put .'s in variable names
pure.data <- read.table("/path/to/pure/data.file")

# a "normal" .tsv with headers and row names
tsv <- read.table("/path/to/tsv/file.tsv", header=T
                  row.names=1, sep="\t")
```

That all worked, so the next question is what is in `tsv` now? A few good methods for investigating `R` objects:
```R
summary(tsv)
length(tsv)
typeof(tsv)
head(tsv)
```

After running these we seem to have a kind of data table known as an `R` Dataframe (which is what read.table returns). To get a better idea, lets explicitly make a data frame, using the `data.frame` command. You may be noticing the `c()` construct popping up, this just defines a vector in R.

```R
test <- data.frame(a=c(1,2,3), b=c("a", "b", "c"))

# In the console, we can print a variable by just typing it
test
  a b
1 1 a
2 2 b
3 3 c

# Dataframes must be rectangular
bad <- data.frame(a=c(1,2,3), b=c("a", "b", "c","1"))

Error in data.frame(a = c(1, 2, 3), b = c("a", "b", "c", "1")) : 
  arguments imply differing number of rows: 3, 4
```

R uses `$` as a sort of operator to reference Dataframe columns by name. We can access each column like so:

```R
test$a
[1] 1 2 3
test$b
[1] a b c
```

We can also use the `[]` notation to references pieces of the frame. This can get quite complicated but in the simpliest way we pass a row and column selector. But first lets add some useful row labels. Often times we will be modifying the row and column names of a Dataframe, so here are a few more useful commands:

```R
# This print the current values
rownames(test)
colnames(test)

# we can also set them
# make sure the dimensions match!
rownames(test) <- c("first", "second", "third")
rownames(test)
       a b
first  1 a
second 2 b
third  3 c
```

Now that we have some nicer row names, lets go back to the `[]` notation. We specify which rows and columns we'd like to "collect" from. Rows first and then columns, seperated by a comma.
```R
# whole thing again
test

# the element "1" which is "first","a"
test["first", "a"]

# what if we want a whole row?
test["first"]
Error in `[.data.frame`(test, "first") : undefined columns selected

# we need to specify the columns we want (all in this case)
test["first",]
      a b
first 1 a

# same goes with a whole column
test[,"a"]
[1] 1 2 3

# We can also use vectors to select from multiple row/cols
test[c("first","third"),]
      a b
first 1 a
third 3 c
```

That is a basic run down of Dataframes. There is a whole lot more too them but we will learn a bit more later. One last thing: we unsuprisingly will want to write our resulting data to a file using the `write.table` function.

```R
?write.table

Usage:

     write.table(x, file = "", append = FALSE, quote = TRUE, sep = " ",
                 eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                 col.names = TRUE, qmethod = c("escape", "double"),
                 fileEncoding = "")
```
