# IBDrecomb Manual 

Content
-------

- [Introduction](#introduction)
- [Citation](#citation)
- [Installation and Usage](#installation-and-usage)
- [License](#license)

# Introduction

**IBD** based **recomb**ination estimation (**IBDrecomb**) is an efficient method for constructing fine-scale genome-wide recombination maps. 
It can be applied to IBD segments inferred from population-based samples with thousands of individuals. IBD segments may be inferred using Refined IBD (add link).

# Citation

TBD

# Installation and Usage

To install IBDrecomb in a linux environment, download the source files and decompress them if necessary, then go to the directory 'src/' and type
```
make
```
the executable file __IBDrecomb__ can be used to estimate recombination rates based on the input IBD segments and customized parameters.
For example, you can go to the example directory and run test as
```
cd ../example/
../src/IBDrecomb -refinedibd test.1.ibd.gz -fbin 10000 > test.10000.map
```
The IBD input can be a plain text file or a gzip-compressed file (.gz), which can be in refined-IBD's format ('-refinedibd') or a generic format ('-ibd') that each line recording the starting position and ending position of an IBD segment. 
In this example, the parameter '-fbin 10000' specifies the estimation at scale 10kb and the output is saved to the file 'test.10000.map'. 
In reality, the estimation interval might be a slight smaller than the customized size because the whole estimation region length (in bp) is not capable of being divided by the interval size without a remainder.
The default output map is in plink map format, which includes four columns: chromosome name, variant identifier (set as '.'), genetic position, and physical coordinate.

Recombination rates and genetic positions are normalized by the parameter '-rate', 
which is the average recombination rate over the estimation region. The value of the rate parameter should be obtained from a family-based reference map.
If the average rate is 1.21cM/Mb, then we can run the program as

```
../src/IBDrecomb -refinedibd test.1.ibd.gz -fbin 10000 -rate 1.21 > test.10000.normalized.map
```
 
Users can also customize advanced parameters such as the bin size. 
A full list of parameters and their default values can be found by typing 

```
../src/IBDrecomb -help

  Usage: IBDrecomb [options] parameters

 (Required inputs:)
 	-refinedibd/-ibd <filename>	# IBD segments input, '-refinedibd' for the refined-IBD's output format and '-ibd' for the two-column generic format
 	-fbin <integer> # Bin size(bp) for outputing recombination rates

 (Optional parameters:)
 	-bin 500000 # Bin size(bp) for first step of estimation
 	-rate 1.0 # Average cM/Mb to normalize the estimated recombination rates
 	-trim 1 # Length(bp) to trim at chromosome ends
 	-fold 1 # Relative size of adjunct-region to the end-region (see Methods of the published paper)
 	-iter 20 # Number of iterations 
	-chr . # Chromosome to estimate
	-help # Print this help file
```

# License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

[\[top\]](#content)
