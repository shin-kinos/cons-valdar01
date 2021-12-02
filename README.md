# cons-valdar01  

A Rust program that scores residue conservation a site in a Multiple Sequence Alignment (MSA) by using Sum-of-pairs measure (SP measure). 

# Description 

* The program scores residue conservation in a MSA. 
* The scoring method is based on SP measure. 
* It takes account of normalized substitution matrices, sequence weighting and gap penalty. 

## Dependencies 

* `colored` ( https://github.com/mackwic/colored ) 

``` 
[dependencies]
colored = "2.0"
``` 

## Installation 

You can compile this program by using `Cargo`. 📦🦀

[e. g.] 

``` 
% cd cons-valdar01-main
% cargo build --release
``` 
Then the object file is generated in `./target/release` directory. 

## Scoring method 

### Conservation score 
The conservation score is calculated based on SP measure as follows [1]: 

<img width="1437" alt="equation_01" src="https://user-images.githubusercontent.com/83740080/142534293-12470ae4-a2c9-4f01-ac8b-783b67fe3aed.png">  

where ***N*** is the length of a site (or number of the sequences in the MSA ), ***Sx( i )*** is an amino acid at site ***i*** in sequence ***x*** and ***M*** is a normalized substitution matrix. 

The normalized substitution ***M*** is given based on Linear method as follows: 

<img width="1440" alt="equation_02" src="https://user-images.githubusercontent.com/83740080/142535042-4ecd251f-3044-4387-8ab4-b053b032afa4.png">   

where ***m*** is a typical substitution scoring matrix such as BLOSUM series. ***min( m )*** and ***max( m )*** are minimum and maximum elements of the substitution scoring matrix respectively. If ***a*** or ***b*** is gap, ***M(a, b)*** gets zero as gap penalty. 

And ***λ*** is for bounding the scores from zero to one: 

<img width="1440" alt="equation_03" src="https://user-images.githubusercontent.com/83740080/142522125-43555ca9-f61c-485b-abd1-efcd5054a8a9.png"> 

Where ***Wx*** is a weighting factor of sequence ***x***. 

### Substitution scoring matrices 

This program supports 11 substitution scoring matrices as follows: 

* BLOSUM45 [2]
* BLOSUM50 
* BLOSUM62 
* BLOSUM70 
* BLOSUM80 
* BLOSUM90 
* PAM30 
* PAM70 
* PAM250 [3] 
* Modified version of PET91* [4]
* Modified version of BLOSUM62*

*These matrices are modified so that all of their diagonal elements are constant at the maximum (`W vs W = 11` in BLOSUM62 and `W vs W = 15` in PET91). 

### Sequence weighting 

This program supports 2 types sequence weighting, the Position-Based method by Henikoff-Henikoff [5] and the Distance-Based method by Vingron-Argos [6]. 

## Input file format 

Aligned Multi-FASTA format. NOTE that nucleotide sequences are not supported. 

See some example input files in `demo` directory. 

## Usage 

Major arguments:

`-i` : Input filename in aligned Multi-FASTA format, REQUIRED.

`-o` : Output filename, REQUIRED. 

`-w` : Methods of sequence weighting ("hen" or "va", default "hen"). 

`-m` : Substitution scoring matrices (default "blosum62"). 

[e. g.]

``` 
% ./cons-valdar01 -i input.fasta -o output.txt -m pet91mod -c yes -t no
``` 

Type `-h` to see other available options. 

## Output 

Site number `\t` Conservation score `\t` Composition of the site

[e.g.] 

<img width="888" alt="readme_result_valdar01" src="https://user-images.githubusercontent.com/83740080/144342434-eb17dc05-b9c5-4ad9-aaf5-1867e27f3d4e.png">　

## References 

1. Valdar, William SJ. "Scoring residue conservation." Proteins: structure, function, and bioinformatics 48.2 (2002): 227-241. 
2. Henikoff, Steven, and Jorja G. Henikoff. "Amino acid substitution matrices from protein blocks." Proceedings of the National Academy of Sciences 89.22 (1992): 10915-10919. 
3. Dayhoff, M., R. Schwartz, and B. Orcutt. "A model of evolutionary change in proteins." Atlas of protein sequence and structure 5 (1978): 345-352. 
4. Jones, David T., William R. Taylor, and Janet M. Thornton. "The rapid generation of mutation data matrices from protein sequences." Bioinformatics 8.3 (1992): 275-282. 
5. Henikoff, Steven, and Jorja G. Henikoff. "Position-based sequence weights." Journal of molecular biology 243.4 (1994): 574-578. 
6. Vingron, Martin, and Patrick Argos. "A fast and sensitive multiple sequence alignment algorithm." Bioinformatics 5.2 (1989): 115-121. 
