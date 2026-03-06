# SeqTrans

Python Tool for multi frame nucleic acid translation and ORF extraction into protein amino acid sequences.

usage: seqTransV4.py [-h] [--output_fasta DOUT] [--3frame] [--no_isobaric] [--longest_orf] [--min_orf_size OSIZE] [--startcodon] [--ncstartcodon] [--ignoregeneid] [--biotype BIOTYPE] [--zipped] *.fasta|*.fa|*.fq|*.fastq|*.gz

Read in na FASTA or FASTQ file and translate to protein sequences. james.wright@icr.ac.uk 2019

positional arguments:
  *.fasta|*.fa|*.fq|*.fastq|*.gz    FASTA or FASTQ file to be translated.

options:
  -h, --help                        show this help message and exit
  --output_fasta DOUT, -o DOUT      Set file to write proteins. Default=TranslatedProteins.fa
  --3frame, -f                      Translated in 3 frames. Default=false
  --no_isobaric, -i                 Do not make proteins isobaric. Default=false
  --longest_orf, -l                 Only store the longest ORF for each Frame. Default=false
  --min_orf_size OSIZE, -m OSIZE    Set minimum amino acid length of ORFs saved. Default=10
  --startcodon, -s                  Require ORFs to have canonical start codon. Default=false
  --ncstartcodon, -n                Require ORFs to have canonical or noncanonical start codon. Default=false
  --ignoregeneid, -g                Do not include gene or transcript id in header. Default=false
  --biotype BIOTYPE, -b BIOTYPE     Add Biotype to Header. Default=Unknown
  --zipped, -z                      Is input file compressed (gzip)? Default=false
