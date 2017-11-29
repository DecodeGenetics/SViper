SViper
=======

This pipeline polishes deletions and isnertions called on long read data (ONT) using short exact reads for refinement.

Installation
------------

For installation, simply clone the repo and use make to compile any utilities.

~~~~
~$ git clone git@lsource2.decode.is:svenjam/SViper.git
~$ cd SViper
~$ make
~~~~

Note: You need a compiler that supports c++14.

Dependencies
------------

Currently, the polishing pipeline lacks a convinient long read pairwise alignment function. Therefore, the main program outputs a fasta file that needs to be mapped by a long read mapper of your choice.

* Linux operating system
* Code so far only tested for gcc 5.4.0
* C++14 support
* A long read mapper of your choice. Recommended: _minimap2_

- - - -

Using SViper
---------------

You can look at all the input requirements by calling the sviper help page

~~~~
~$ sviper -h
~~~~

The usage of the whole polishing pipeline, needs a remapping step with a long read mapper of your choice, though.
An example of using sviper is the following:

#. Call sviper
~~~~
~$ sviper -s short-reads.bam -l long-reads.bam -r ref.fa -c variants.vcf -o example -g example.log
~~~~
This will output a `example.fa` file, that contains all the polished sequences for recalling refined variants.

#. Map polished sequences to reference
~~~~
~$ my-fav-mapper -input example.fa -r ref.fa
~~~~
This will output a sam or bam file that serves as input for the third step.

#. Evaluate the mapping and create a **polished VCF file**
~~~~
~$ evaluate_final_mapping my-fav-mapper-output.bam variants.vcf
~~~~
This will output a file called `variants.vcf.polished.vcf` containing variants with the same ID but refined break points or even failures if the variant was "polished away" (see Interpreting the output).

### IMPORTANT

There are several requirements for using the polishing:

#. The vcf file must be a structural variant format (tags instead of sequences, e.g. <DEL>). ALso the INFO field must include the END tag, giving the end position of the variant, as well as the SVLEN tag in case of insertions.
#. The bam files must be indexed.
#. The reference sequence (FASTA) must be indexed.
#. (Obvisoulsy, the bam files should correspond to the same individual/sample mapped to the given reference.)

### Input arguments:

* `-c, --candidate-vcf` The candidate vcf file
    This file contains all the variants that will be polished in the run.
    If you want to accelarate polishing on a large vcf file, you might want to split the file into several small ones and call sviper each one seperately.

* `-s, --short-read-bam` (BAM_FILE)
          The indexed bam file containing short used for polishing at variant sites.

* `-l, --long-read-bam` (BAM_FILE)
          The indexed bam file containing long reads to be polished at variant sites.

* `-r, --reference` (FA_FILE)
          The indexed (fai) reference file.

* `-o, --output-fa` (FA_FILE)
          A filename for the output file. NOTE: The current output is a fasta file, that contains the polished sequences for each variant. Since the final realignment is not part of this tool yet,The user must map the fasta file with a mapper of his choice (e.g. minimap2) and then call evaluate_final_alignment.

* `-g, --log-file` (TXT_FILE)
          A filename for the log file. Default: polishing.log.

* `-k, --flanking-region` (INT)
          The flanking region in bp's around a breakpoint to be considered for polishing In range [50..1000]. Default: 400.

            ~~~~
                            start x             end y
            ------------------|------------------|----------------
                       !------------!      !------------!
                      x-50        x+50    y-50         y+50
            ~~~~

### Contact
If you have any questions write me a mail: svenja.mehringer[AT]gmail.com
