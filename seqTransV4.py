#!/usr/bin/python2
# Sequence Translator (Foward Frames Only)
# James Wright 2019

#used to get cmd line arguments 
import argparse
import re
import gzip
import glob


#Read command line arguments and create help documentation using argparse
parser = argparse.ArgumentParser(
	description='''Read in na FASTA or FASTQ file and translate to protein sequences. james.wright@icr.ac.uk 2019''')
parser.add_argument('fasta', metavar='*.fasta|*.fa|*.fq|*.fastq|*.gz', help='FASTA or FASTQ file to be translated.')
parser.add_argument('--output_fasta', '-o', dest='dout', default='TranslatedProteins.fa', help='Set file to write proteins. Default=TranslatedProteins.fa')
parser.add_argument('--3frame', '-f', dest='frames', default=False, action='store_true', help='Translated in 3 frames. Default=false')
parser.add_argument('--no_isobaric', '-i', dest='iso', default=False, action='store_true', help='Do not make proteins isobaric. Default=false')
parser.add_argument('--longest_orf', '-l', dest='longest', default=False, action='store_true', help='Only store the longest ORF for each Frame. Default=false')
parser.add_argument('--min_orf_size', '-m', type=int, dest='osize', default=10, help='Set minimum amino acid length of ORFs saved. Default=10')
parser.add_argument('--startcodon', '-s', dest='scodon', default=False, action='store_true', help='Require ORFs to have canonical start codon. Default=false')
parser.add_argument('--ncstartcodon', '-n', dest='ncscodon', default=False, action='store_true', help='Require ORFs to have canonical or noncanonical start codon. Default=false')
parser.add_argument('--ignoregeneid', '-g', dest='ignoregene', default=False, action='store_true', help='Do not include gene or transcript id in header. Default=false')
parser.add_argument('--biotype', '-b', dest='biotype', default='Unknown', help='Add Biotype to Header. Default=Unknown')
parser.add_argument('--zipped', '-z', dest='uzip', default=False, action='store_true', help='Is input file compressed (gzip)? Default=false')

#parser.add_argument('--id_filter', '-x',  dest='flist', default='', help='A file containing a list of IDs to include in the results. Default=NA (include all ids)')

args = parser.parse_args()


#Basic genetic code to use for translation
gcode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
      'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

ncscodons = ['CTG','TTG','GTG','ACG','ATA','ATT','ATC']
stopcoons = ['TAA','TAG','TGA']

#Function to translate na into aa sequence
def translate ( se ):
      translate = ''.join([gcode.get(se[3*i:3*i+3],'X') for i in range(len(se)//3)])
      return translate

def orf_translate ( se, sepos, scodon, ncscodon, iso):
	
	#Init dictionary of translated sequences and positions
	torfs = {}

	#Intitialize empty orf tranlated sequence
	current_torf = ''

	#Set starting position
	spos = sepos

	#Use flags to track if start codon has been found if required
	start_found = True
	if scodon == True or ncscodon == True:
		start_found = False

	#Loop each codon in sequence
	for pos in range(0, len(se)-2, 3):
		codon = se[pos:pos+3]

		#Check for canonical start codon if required and update start position
		if start_found == False and codon == 'ATG':
			start_found = True
			spos = sepos + pos

		#Check for non-canonical start codons if required and update start position	
		if start_found == False and ncscodon == True and codon in ncscodons:
			start_found = True
			spos = sepos + pos
			
		#Check for stop codon
		if codon in stopcoons:

			#If start codon found store current orf in dictionary with starting position
			if start_found == True and current_torf != '':

				#Make isobaric if selected
				if iso == True:
					current_torf = current_torf.replace('I','L')

				torfs[spos] = current_torf

			#Reset current orf sequence and update start position
			current_torf = ''
			spos = sepos + pos + 3
			if scodon == True or ncscodon == True:
				start_found = False

		else:

			#If start codon found append aa to current orf sequence
			if start_found == True:
				current_torf += gcode.get(codon,'X')

	return torfs



#Function to process each FASTA sequence by selected frame(s) and orf filtering parameters
def process_sequence ( seq, id, nframes, iso, osize, longest, oFa, scodon=False, ncscodon=False, ingnoregene=False, biotype='Unknown' ):

	olen = 0	#longest ORF length
	lf = 99		#longest ORF frame
	lc = -1		#longest ORF id
	lorf =''	#longest ORF sequence
	lspos = -1 	#longest ORF transcript start position
	lcnt = 1	#counter for ORF enumeration

	#split id on any white space so acession can be modified
	idsp = id.split(None, 1)

	#Extract gene ID by removing '>' from start of FASTA header
	grid = idsp[0].replace('>','')

	#Extract transcript ID by removing spladder isoform info if present 
	trid = re.sub('_iso\d+', '', grid)

	#Set gene ID for FASTA header
	gid = grid
	if len(idsp) > 1: gid = idsp[1]

	#Loop each forawrd frame (1 or all 3)
	for frame_start in range(nframes):

		#Get sequence for this frame
		fseq = seq[frame_start:]
	
		#Get dictionary of translated ORFs and their start positions
		orfs = orf_translate(fseq, frame_start, scodon, ncscodon, iso)	

		#Loop and enumerate each ORF found get pair of position and orf sequence
		for ocnt, (start_pos, orf) in enumerate(orfs.items()):

			#Does the orf meet minimum length requirement from args
			if len(orf) >= osize:

				#if option to only retain longest ORF is set check current orf against longest so far and update stored information
				if longest == True:

					if (len(orf) > olen):
						olen = len(orf)
						lorf = orf
						lf = frame_start
						lc = ocnt
						lspos = start_pos
						
				#>alt_3prime.4188_iso1|1 frame=0 offset=0 orf=0 gene=ENSG00000188976.9 geneID=alt_3prime.4188 transcriptID=alt_3prime.4188_iso1
				#>alt_3prime.4188_iso1_Frame0_spos0 orf=0 gene=ENSG00000188976.9

				#Otherwise write each translated ORF into new FASTA
				else:
					oFa.write(idsp[0] + '|' + str(lcnt) + ' frame=' + str(frame_start) + ' offset=' + str(start_pos) + ' orf=' + str(ocnt))
					if ingnoregene == False:
						oFa.write(' gene=' + gid  + ' geneID=' + grid + ' transcriptID=' + trid )
					if biotype != '' or biotype != 'Unknown':
						oFa.write(' biotype=' + biotype)
					oFa.write('\n')
					oFa.write(orf + '\n')
					lcnt+=1

	#If selected write longets ORF to FASTA after looping all ORFs in all Frames
	if longest == True:
		#Does the longest ORF meet minimum length requirements
		if len(lorf) >= osize:
			oFa.write(idsp[0] + '|' + str(lcnt) + ' frame=' + str(lf) + ' offset=' + str(lspos) + ' orf=' + str(lc) )
			if ingnoregene == False:
				oFa.write(' gene=' + gid  + ' geneID=' + grid + ' transcriptID=' + trid )
			if biotype != '' or biotype != 'Unknown':
				oFa.write(' biotype=' + biotype)
			oFa.write('\n')
			oFa.write(lorf + '\n')

	return

#Open zipped or unzipped fastq files and return handle
def open_fasta (file, zipped):
	if zipped:
		print ("ZIPPED")
		return (gzip.open(file, 'r'))
	else:
		print ("NORMAL")
		return( open(file, 'r'))

#################################################

#empty protein sequence, protein id
sequ = ''	
proID = ''

#Set number of frames to be translated
nf = 1
if args.frames == True:
	nf = 3

#print (parser.parse_args())	

#Open FASTA file for writing translated sequences
outFasta = open(args.dout, 'w')

#Open FASTA file using first cmd line argument
fasta = open_fasta(args.fasta, args.uzip)

#loop each line in the file
for line in fasta:

	#Uncompress line if zipped
	if args.uzip:
		line = line.decode()

	#if this line starts with ">" then process sequence if not empty
	if line[0] == '>':
		if sequ != '':
			#Process and translate current sequence
			process_sequence(sequ.upper(), proID, nf, args.iso, args.osize, args.longest, outFasta, args.scodon, args.ncscodon, args.ignoregene, args.biotype)

		#Extract next sequence identifier (FASTA Header)
		proID = line.rstrip()
		#New sequence
		sequ = ''

	elif line[0] == '@':	#if this is a fastq file and line starts with @ then extract sequence identifier
		if sequ != '':
			#Process and translate current sequence
			process_sequence(sequ.upper(), proID, nf, args.iso, args.osize, args.longest, outFasta, args.scodon, args.ncscodon, args.ignoregene, args.biotype)

		#Extract next sequence identifier (FASTQ Header)
		proID = (line.rstrip().split())[0].replace('@','')
		#New sequence
		sequ = ''

	#if not accession line and not quality line then append aa sequence (with no newline or white space) to seq string
	elif line.rstrip().isalpha(): 
		sequ+=line.rstrip()
		
#Close files
fasta.close()

#At end of FASTA file process final sequence
if sequ != '':
	process_sequence(sequ.upper(), proID, nf, args.iso, args.osize, args.longest, outFasta, args.scodon, args.ncscodon, args.ignoregene, args.biotype)


outFasta.close()
