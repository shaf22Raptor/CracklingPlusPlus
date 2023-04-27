# CracklingPlusPlus

Rapid Whole-Genome Identification of High Quality CRISPR Guide RNAs with the Crackling Method


## Preamble

> The design of CRISPR-Cas9 guide RNAs is not trivial. In particular, evaluating the risk of off-target 
modifications is computationally expensive: a brute-force approach would require comparing each candidate 
guide with every possible CRISPR target site in the genome. In a mammalian genome, this means hundreds of 
millions of comparisons for each guide. We have previously introduced Crackling, a gRNA design tool that 
relies on Inverted Signature Slice Lists (ISSL) to accelerate off-target scoring by only considering sites
with partial matches (a slice) with the candidate guide. This produced an order of magnitude speed up whilst 
still maintaining scoring accuracy. Here, we present a complete reimplementation of Crackling in C++ and 
discuss further improvements. Using longer slices means fewer comparisons, and we show it is possible to 
construct a collection of slices that still preserve an exact off-target score. We have benchmarked two ISSL
configurations with the new version of Crackling and report a 15-22 times speed up over the default ISSL 
configuration. This increased performance comes at the cost of increased memory usage, but this can be offset
by using memory mapped files, and we show that this has no significant impact on performance. 
CracklingPlusPlus is available at https://github.com/bmds-lab/CracklingPlusPlus under the Berkeley Software 
Distribution (BSD) 3-Clause license.

## Dependencies
- [CMake](https://cmake.org/)
- [Boost](https://www.boost.org/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [RNAfold](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)

Note: Please ensure that your version of Boost that you have installed is supported by your version of CMake.

## Installation
1. Clone or download the repo.
```bash
git clone https://github.com/bmds-lab/CracklingPlusPlus.git ~/CracklingPlusPlus
cd ~/CracklingPlusPlus
```

2. Create build directory
```bash
mkdir build
cd build
```

3. Run CMake to generate build files

```bash
CMake ..
```

4. Run build system command. E.g. `make`
```bash
make
```
All of the programs (CracklingPlusPlus, ISSLCreateIndex and ExtractOfftargets) have now been built.

## Building Bowtie2 Index
The Bowtie2 manual can be found [here](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

Our recommended usage:
```bash
bowtie2-build --threads 128 input-file output-file
```

For example:
```bash
bowtie2-build --threads 128 ~/genomes/mouse.fa ~/genomes/mouse.fa.bowtie2
```

Note: Bowtie2 produces multiple files for its index. When setting the Bowtie2 index variable in the `config.ini` file please use the value you used for `output-file`. So for the example above, you would set `bowtie2-index = ~/genomes/mouse.fa.bowtie2`
## Extract off-targets for ISSL Index
Note: You will need to ensure you have completed the installation step before completing this step as the installation will also install the program needed to extract off-targets. 

The ExtractOfftargets binary can be found in build folder. Based on the install instructions that will be:
```
~/CracklingPlusPlus/build/ExtractOfftargets/ExtractOfftargets
```

Usage:
```bash
ExtractOfftargets <output-file>  {<input-file-1> <input-file-2> ... <input-file-n> | <input-dir>}
```

Arguements:
```
output-file: A filepath to save the off-target sites

input-file-x: A single, or a space serpated list, of multi-FASTA formatted files

input-dir: A directory, containing multiple multi-FASTA formatted files. (Note: This will process EVERY file in the directory)
```

For example:
```bash
ExtractOfftargets ~/genomes/mouse_offtargets.txt ~/genomes/mouse.fa
```
or
```bash
ExtractOfftargets ~/genomes/mouse_offtargets.txt ~/genomes/mouse_chromosone_1.fa ~/genomes/mouse_chromosone_2.fa ~/genomes/mouse_chromosone_3.fa
```
or
```bash
ExtractOfftargets ~/genomes/mouse_offtargets.txt ~/genomes/mouse_fasta_files/
```

## Building ISSL Index

Note: You will need to ensure you have completed the installation step before completing this step as the installation will also install the program needed to build the ISSL Index. 

The ISSLCreateIndex binary can be found in build folder. Based on the install instructions that will be:
```
~/CracklingPlusPlus/build/ISSLCreateIndex/ISSLCreateIndex
```

Usage:
```bash
ISSLCreateIndex <offtarget-sites> <slice-config> <sequence-length> <output-file>
```
Arguements:
```
offtarget-sites: A text file containing off-target sites

slice-config: A text file containing a set of slice configurations (See samples folder in repository)

sequence-length: The length of an off-target site

output-file: A filepath to save the ISSL index
```

For example:
```bash
ISSLCreateIndex ~/genomes/mouse_offtargets.txt ~/CracklingPlusPlus/sample/slice4-5.txt 20 ~/genomes/mouse_indexed.issl
```
## Running CracklingPlusPlus
Please ensure all of the above steps have been completed before running the program. 
To run the program simply fill out the provided `config.ini` in the samples folder and call the program as follows:

```bash
CracklingPlusPlus <config-file>
```



## References

Ben Langmead and Steven L Salzberg. Fast gapped-read alignment with Bowtie2. Nature Methods, 9(4):357, 2012.

Bradford, J., Chappel, T., & Perrin, D. (2022). Rapid Whole-Genome Identification of High Quality CRISPR Guide RNAs with the Crackling Method. The CRISPR Journal, 5(3), 410-421.

Bradford, J., & Perrin, D. (2019). A benchmark of computational CRISPR-Cas9 guide design methods. PLoS computational biology, 15(8), e1007274.

Bradford, J., & Perrin, D. (2019). Improving CRISPR guide design with consensus approaches. BMC genomics, 20(9), 931.

Chari, R., Yeo, N. C., Chavez, A., & Church, G. M. (2017). sgRNA Scorer 2.0: a species-independent model to predict CRISPR/Cas9 activity. ACS synthetic biology, 6(5), 902-904.

Lorenz, R., Bernhart, S. H., Zu Siederdissen, C. H., Tafer, H., Flamm, C., Stadler, P. F., & Hofacker, I. L. (2011). ViennaRNA Package 2.0. Algorithms for molecular biology, 6(1), 1-14.

Montague, T. G., Cruz, J. M., Gagnon, J. A., Church, G. M., & Valen, E. (2014). CHOPCHOP: a CRISPR/Cas9 and TALEN web tool for genome editing. Nucleic acids research, 42(W1), W401-W407.

Sunagawa, G. A., Sumiyama, K., Ukai-Tadenuma, M., Perrin, D., Fujishima, H., Ukai, H., ... & Shimizu, Y. (2016). Mammalian reverse genetics without crossing reveals Nr3a as a short-sleeper gene. Cell reports, 14(3), 662-677.