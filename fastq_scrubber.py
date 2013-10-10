#!/usr/bin/python

###########################
# E. Roberson             #
# Created 2013-October-07 #
###########################

###########
# Imports #
###########
from datetime import datetime
from collections import namedtuple
import re
import sys

####################
# Version and name #
####################
SCRIPT_PATH = sys.argv[0]
SCRIPT_NAME = SCRIPT_PATH.split('/|\\')[-1]
VERSION = '1.0.0'

#######################
# Constants & globals #
#######################
CleanResult = namedtuple( 'CleanResult', 'sequence quality goodBp incr_adapt good_seq' )
SequenceInterpretation = namedtuple( 'SequenceInterpretation', 'totalBp adjustBp extensionName' )

# default adapters
DEFAULT_ADAPTERS = ('AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC', 'GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA', 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG', 'GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG')
ADAPT_SEQ = ''

# Misc constants and counts
MIN_SEQ_SIZE = 5
N_PERC_CUTOFF = 0.85
MIN_CMD_ARGS = 3

# Quality score stuff
ILMN_OFFSET = 64
SNGR_OFFSET = 33

illumina2phredDict = {}
for i in range(-5,42):
	illumina2phredDict[chr(i+ILMN_OFFSET)] = chr(i+SNGR_OFFSET)

phred2illuminaDict = {}
for i in range(0,42):
	phred2illuminaDict[chr(i+SNGR_OFFSET)] = chr(i+ILMN_OFFSET)

# Alignment to adapter
MATCH_SCORE = 1
MISMATCH_SCORE = 0
SUBSTITUTION_MATRIX = {'A':{'A':MATCH_SCORE, 'G':MISMATCH_SCORE,'C':MISMATCH_SCORE, 'T':MISMATCH_SCORE, 'N':MATCH_SCORE}, 'G':{'A':MISMATCH_SCORE, 'G':MATCH_SCORE,'C':MISMATCH_SCORE, 'T':MISMATCH_SCORE, 'N':MATCH_SCORE}, 'C':{'A':MISMATCH_SCORE, 'G':MISMATCH_SCORE,'C':MATCH_SCORE, 'T':MISMATCH_SCORE, 'N':MATCH_SCORE}, 'T':{'A':MISMATCH_SCORE, 'G':MISMATCH_SCORE,'C':MISMATCH_SCORE, 'T':MATCH_SCORE, 'N':MATCH_SCORE}, 'N':{'A':MATCH_SCORE, 'G':MATCH_SCORE,'C':MATCH_SCORE, 'T':MATCH_SCORE, 'N':MATCH_SCORE}}

############
# Counters #
############
readCount = 0
badReadCount = 0
highNsCount = 0
escapesOfDeath = 0
adapterTrimmedCount = 0
totalSequenceCount = 0
usableSequenceCount = 0

###########
# Globals #
###########
QUALITYSCALE = ""
UNITTEST = False
RIGHTQUALTRIM = True
junkHandles = []
goodHandles = []

####################################
# Early fxn definition for classes #
####################################
def nonZero( count ):
	return count if count > 0 else 1

def errMsg( s ):
	print s
	sys.exit(1)

#####################
# class Definitions #
#####################
class SearchAdapter:
	def __init__(self, seq, start, length):
		self.subSeq = seq[start:(start+length)]
		self.start = start
		self.size = length
		self.sequence = seq
		self.regex = re.compile(self.subSeq)
		self.sequenceSize = len(seq)
		
	def adapterSearch( self, query ):
		return self.regex.search( query )

###################
# Test cline args #
###################
usage = """
%s v%s
Usage: python %s <infile> <outfiles> <optional_arguments>

REQUIRED
========
<infile>      Name (current directory) or path to raw read file(s)
              e.g. s_x_1.fastq,s_x_2_.fastq
<outfile>     Output file base name
              Outname will be outname_1|2_clean.fq|fqs & outname_1|2_junk.fq|fqs

OPTIONAL
========
outscale=val  Use 'phred', 'illumina1.3' or 'illumina1.0'
r1_ltrim=#    Read 1, hardclip # bases from left end
r1_rtrim=#    Read 1, hardclip # bases from right end
r2_ltrim=#    Read 2, hardclip # bases from left end
r2_rtrim=#    Read 2, hardclip # baes from the right end
stripAdapter  Strips adapter from 3' read ends*
adapter=N..N  Comma-separated adapter sequence to trim from 3' ends
noRightQtrim  Don't trim off low quality reads from 3' end
trimCutoff=#  Trim bases off 3' end with quality lower than # (not offset)

*if no adapters are specified with the 'adapter=NN,NN' command
 the default adapters in this file are stripped.

Example:
python %s s_6_1.fq,s_6_2.fq s_6 r1_rtrim=3 stripAdapter noRightQtrim
""" % (SCRIPT_PATH, VERSION, SCRIPT_NAME, SCRIPT_NAME)

##################################
# check for valid command number #
##################################
argc = len(sys.argv)

if 'unittest' in sys.argv:
	UNITTEST = True

else:
	if argc < MIN_CMD_ARGS:
		errMsg( "Too few arguments. Specify appropriately or select \'unittest\' to test script functionality\n%s\n" % (usage) )

##################################
# Set up important global values #
# and defaults                   #
##################################
if not UNITTEST:
	INFILE = sys.argv[1].split(',')
	OUTFILE = sys.argv[2]

PHRED_MIN = 33
ILMN1P3_MIN = 64
ILMN1P0_MIN = 59
OUTPUTEXTENSION = '.fqs'
JNKEXT = '.scarf'
OUTPUTSCALE = 'phred'
LOW_SCORE = chr( PHRED_MIN )
R1LTRIM = 0
R1RTRIM = 0
R2LTRIM = 0
R2RTRIM = 0
REM_ADAPT = False
RE_LEN = 6
RE_NUM = 2
RE_STEP = 4
ADAPT_SCORE_CUTOFF = RE_LEN * MATCH_SCORE
SEARCH_LIST = []
TRIMCUTOFF = 5
ENDS = 'paired'
INTYPE = 'fastq'

############################
# Parse optional arguments #
############################
if argc > MIN_CMD_ARGS:
	for i in range(MIN_CMD_ARGS, argc):
		try:
			(argkey, argvalue) = sys.argv[i].split('=')
		except:
			argkey = sys.argv[i]
		
		if argkey == 'outscale':
			if argvalue.lower() == "phred":
				OUTPUTSCALE = 'phred'
				OUTPUTEXTENSION = ".fqs"
				LOW_SCORE = chr( PHRED_MIN )
			elif argvalue.lower() == "illumina1.3":
				OUTPUTSCALE = 'illumina1.3'
				OUTPUTEXTENSION = ".fq"
				LOW_SCORE = chr( ILMN1P3_MIN )
			elif argvalue.lower() == "illumina1.0":
				OUTPUTSCALE = 'illumina.0'
				OUTPUTEXTENSION = ".fqlod"
				LOW_SCORE = chr( ILMN1P0_MIN )
			else:
				errMsg( "scale value [%s] not recognized" % (argvalue) )
		elif argkey == "r1_ltrim":
			try:
				R1LTRIM = int(argvalue)
			except:
				errMsg( "Read 1, left trim value of %s not recognized as an integer" % (argvalue) )
		elif argkey == "r1_rtrim":
			try:
				R1RTRIM = int(argvalue)
			except:
				errMsg( "Read 1, right trim value of %s not recognized as an integer" % (argvalue) )
		elif argkey == "r2_ltrim":
			try:
				R2LTRIM = int(argvalue)
			except:
				errMsg( "Read 2, left trim value of %s not recognized as an integer" % (argvalue) )
		elif argkey == "r2_rtrim":
			try:
				R2RTRIM = int(argvalue)
			except:
				errMsg( "Read 2, right trim avlue of %s not recognized as an integer" % (argvalue) )
		elif argkey == "adapter":
			REM_ADAPT = True
			ADAPT_SEQ = argvalue.upper().split(',')
		elif argkey	== "stripAdapter":
			REM_ADAPT = True
			if len( ADAPT_SEQ ) == 0:
				ADAPT_SEQ = DEFAULT_ADAPTERS
		elif argkey == "noRightQtrim":
			RIGHTQUALTRIM = False
		elif argkey == "qualTrimCutoff":
			TRIMCUTOFF = int( argvalue )
		else:
			errMsg( "Option [%s] not recognized. See usage for appropriate commands" % (argkey) )

########
# Fxns #
########
def setQualScale(filename, inputType=INTYPE if not UNITTEST else None):
	try:
		INFH = open(filename, 'r')
	except:
		errMsg( "Error when trying to open file [%s] to determine quality scale" % (filename) )
	
	# phred		 0 to 93, ASCII 33 to 126; offset 33; phred scale
	# sol1.0	-5 to 62, ASCII 59 to 126; offset 64; log odds
	# sol1.3+	 0 to 62, ASCII 64 to 126; offset 64; phred scale
	# Qsanger = -10 * log10(p)
	# Qsol1.0 = -10 * log10( p / 1-p )
	# 33 - 58 == phred
	# 59 - 63 == sol1.0
	# sol1.3 if no evidence for either
	
	maxReads = 5000
	readCount = 0
	minQual = 10000
	localIlmn1p0Min = 59
	localIlmn1p3Min = 64
	localPhredMin = 33
	
	if inputType == "scarf":
		for line in INFH:
			readCount += 1
			if readCount > maxReads:
				break
				
			line = line.rstrip("\r\n")
			lineData = line.split(":")
			
			if len(lineData[-2]) != len(lineData[-1]):
				# pesky metacharacters...
				errMsg( "We have a problem! Setting qual score but qual and seq have different lengths!!!" )
			
			for qScore in lineData[-1]:
				if ord(qScore) < minQual:
					minQual = ord(qScore)
					
				if ord(qScore) < localIlmn1p0Min:
					readCounts = maxReads
					break
	elif inputType == 'fastq':
		seqBuffer = ''
		qualBuffer = ''
		junkBuffer = ''
		
		while readCount < maxReads:
			junkBuffer = INFH.readline() # @ identifier
			seqBuffer = INFH.readline().rstrip("\r\n") # seq
			junkBuffer = INFH.readline() # + identifier
			qualBuffer = INFH.readline().rstrip("\r\n") # quality
			
			if len(seqBuffer) != len(qualBuffer):
				continue
			
			readCount += 1
			
			for qScore in qualBuffer:
				if ord(qScore) < minQual:
					minQual = ord(qScore)
					
				if ord(qScore) < localIlmn1p0Min:
					readCounts = maxReads
					break
	
	INFH.close()
	
	if minQual < localIlmn1p0Min:
		return "phred"
	elif minQual < localIlmn1p3Min:
		return "illumina1.0"
	else:
		return "illumina1.3"

def qualityConverter(qualities, inscale, outscale, localIlmn2Phred=illumina2phredDict, localPhred2IlmnDict=phred2illuminaDict):
	if inscale == outscale:
		return qualities
	elif inscale == 'illumina1.3' and outscale == 'phred':
		return ''.join([localIlmn2Phred[val] for val in qualities])
	elif inscale == 'phred' and outscale == 'illumina1.3':
		return ''.join([localPhred2IlmnDict[val] for val in qualities])
	else:
		errMsg( "Error in converting scales. Input scale is %s and output scale is %s" % (inscale, outscale) )
	
def stringLeftSub( string, subValue, num ):
	return "%s%s" % (subValue*num, string[num:])
		
def hardClip( seq, qual, side, length):
	if side == "right":
		return (seq[:-length], qual[:-length])
	elif side == "left":
		return (seq[length:], qual[length:])
	else:
		errMsg( "\"side\" command [%s] not recognized. Probably code error." % (side) )

def nRstrip( seq, qual ):
	newSeq = seq.rstrip("N")
	return( newSeq, qual[:len(newSeq)])

def nLstrip( seq, qual, lowQualScore ):
	leftN = len(seq) - len(seq.lstrip('N'))
	if leftN != 0:
		return stringLeftSub( qual, lowQualScore, leftN )
	else:
		return qual

if RIGHTQUALTRIM:
	def trimLowQualBases(seq, qual, qcutoff):
		if len(seq) > 0:
			lowQualCount = 0
			for curr_char in reversed(qual):
				if ord(curr_char) <= qcutoff:
					lowQualCount += 1
				else:
					break
			
			if lowQualCount > 0:
				return (seq[:-lowQualCount], qual[:-lowQualCount])
		return (seq, qual)
else:
	def trimLowQualBases( seq, qual, qcutoff ):
		return( seq, qual )

def eatAdapters(seq, qual, reList, minAdaptMatchLength=RE_LEN, subMatrix=SUBSTITUTION_MATRIX, adaptScoreCutoff=ADAPT_SCORE_CUTOFF, minAdptIdentity=0.85 ):
	if len(seq) >= minAdaptMatchLength:
		for current_re in reList:
			result = current_re.adapterSearch( seq )
			if result:
				adapterMatchStart = result.start() - current_re.start # remember some of these matches are offset. adjust for it before returning substring for matching
				
				seqMatchSubstring = seq[adapterMatchStart:]
				maxPosition = len(seqMatchSubstring) if len(seqMatchSubstring) < current_re.sequenceSize else current_re.sequenceSize
				
				identityScore = 0
				for alignmentIndex in range(maxPosition):
					identityScore += subMatrix[current_re.sequence[alignmentIndex]][seqMatchSubstring[alignmentIndex]]
				
				if identityScore >= int(minAdptIdentity * maxPosition):
					return (seq[:adapterMatchStart], qual[:adapterMatchStart], 1)
	return (seq, qual, 0)

if not UNITTEST:
	def cleanSequence( seq, qual, read, trimVal, lowQualScore, r1ltrim=R1LTRIM, r1rtrim=R1RTRIM, r2ltrim=R2LTRIM, r2rtrim=R2RTRIM, ADAPTER=ADAPT_SEQ, adaptReList=SEARCH_LIST, minSeqSize=MIN_SEQ_SIZE, adapterRemoval=REM_ADAPT ):
		##########################################################
		# NOTE                                                   #
		# Indexes have already been removed from sequences here. #
		##########################################################
		if read == "read1":
			# hard clip
			if r1rtrim > 0:
				seq, qual = hardClip(seq, qual, "right", r1rtrim)
			if r1ltrim > 0:
				seq, qual = hardClip(seq, qual, "left", r1ltrim)
		elif read == "read2":
			# hard clip
			if r2rtrim > 0:
				seq, qual = hardClip(seq, qual, "right", r2rtrim)
			if r1ltrim > 0:
				seq, qual = hardClip(seq, qual, "left", r2ltrim)
		else:
			errMsg( "In fxn cleanSequence the read value of [%s] is not recognized" % (read) )
		
		#############
		# Remove Ns #
		#############
		# right
		seq, qual = nRstrip( seq, qual )
		
		# left
		qual = nLstrip( seq, qual, lowQualScore )
		
		############################
		# Remove low quality bases #
		############################
		seq, qual = trimLowQualBases( seq, qual, qcutoff=trimVal )
		
		###########################################
		# Strip any remaining adapter from 3' end #
		###########################################
		incr_adapt = False
		if adapterRemoval:
			seq, qual, incr_adapt = eatAdapters( seq, qual, reList=adaptReList )
		
		goodSeq = True if len(seq) > minSeqSize else False
		
		goodBp = len(seq) - seq.count('N')
		
		if len(seq) < minSeqSize:
			sizeDiff = minSeqSize - len(seq)
			seq = "%s%s" % (seq, 'N' * sizeDiff)
			qual = "%s%s" % (qual, lowQualScore * sizeDiff)
		
		return CleanResult( seq, qual, goodBp, incr_adapt, goodSeq )

def boolString( boolVal ):
	if boolVal:
		return "True"
	return "False"

def isIntString( s ):
	try:
		int(s)
	except ValueError:
		return False
	return True

def unittest():
	print "%s v%s\n" % (SCRIPT_NAME, VERSION)
	print "Running function unittests"
	
	# randomly generated sequences for testing functions
	TEST_SEQ1 = 'ACGCGATATGAAGCCTGTGACGGTTAGCTGTCGGAGGAAGTTTACCCTAG'
	#TEST_SEQ2 = 'ATTCGAGTGAAGGCTGGATGGTTCATATAGCTCTTCTGAACAAGGCGTCT'
	#TEST_SEQ3 = 'GATTTTCGTTAATACCCGATGTGTGGGACTGACTGCCGTACTGGGCCATG'
	
	QUAL1 = 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWOvJJ]tiJ'
	#QUAL2 = 'lHh|DjKoDxZutPdoTTCxNpnKgHKPlkERCT}zK_yVOwBxVVTMwz'
	#QUAL3 = 'D{Rxv_nms]SqoFRXUdHADUfzziLW^`txB{[bI_GLF`LEsF^O]Y'
	
	testsPassed = 0
	testsFailed = 0
	currStatus = ""
	
	PASS = "Passed"
	FAIL = "Failed"
	
	######################
	# Quality Conversion #
	######################
	print "\nTesting qualityConverter"
	
	# 64-70 [@ABCDEF]. phred>illumina 95-101 [_`abcde].
	
	if qualityConverter( '@ABCDEF', 'phred', 'phred' ) == '@ABCDEF':
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
		
	print "%s qualityConverter (phred>phred) 1/4" % (currStatus)
	
	if qualityConverter( '@ABCDEF', 'illumina1.3', 'illumina1.3' ) == '@ABCDEF':
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s qualityConverter (illumina1.3>illumina1.3) 2/4" % (currStatus)
	
	if qualityConverter( '@ABCDEF', 'phred', 'illumina1.3' ) == '_`abcde':
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s qualityConverter (phred>illumina1.3) 3/4" % (currStatus)
	
	if qualityConverter( '_`abcde', 'illumina1.3', 'phred' ) == '@ABCDEF':
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s qualityConverter (illumina1.3>phred) 4/4" % (currStatus)
	
	#############################################
	# Substituting left side of quality strings #
	#############################################
	print "\nTesting stringLeftSub"
	
	if stringLeftSub( 'ACGCGATATGAAGCCTGTGACGGTTAGCTGTCGGAGGAAGTTTACCCTAG', 'N', 10 ) == 'NNNNNNNNNNAAGCCTGTGACGGTTAGCTGTCGGAGGAAGTTTACCCTAG':
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s stringLeftSub 1/1" % (currStatus)
	
	###############################
	# Test hard clipping sequence #
	###############################
	print "\nTesting hardClip"
	
	tmpSeq, tmpQual = hardClip( TEST_SEQ1, QUAL1, 'right', 5 )
	
	if tmpSeq == 'ACGCGATATGAAGCCTGTGACGGTTAGCTGTCGGAGGAAGTTTAC' and tmpQual == 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWOvJ':
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s hardClip (right) 1/2" % (currStatus)
	
	tmpSeq, tmpQual = hardClip( TEST_SEQ1, QUAL1, 'left', 5 )
	
	if tmpSeq == 'ATATGAAGCCTGTGACGGTTAGCTGTCGGAGGAAGTTTACCCTAG' and tmpQual == 'YKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWOvJJ]tiJ':
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s hardClip (left) 2/2" % (currStatus)
	
	###############
	# strip 3' Ns #
	###############
	print "\nTesting nRstrip"
	
	tmpSeq, tmpQual = nRstrip( 'ACGCGATATGAAGCCTGTGACGGTTAGCTGTCGGAGNNNGNNNNNNNNNN', 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWOvJJ]tiJ' )
	
	if tmpSeq == 'ACGCGATATGAAGCCTGTGACGGTTAGCTGTCGGAGNNNG' and tmpQual == 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXh':
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s nRstrip 1/1" % (currStatus)
	
	###################################
	# strip 5' Ns by quality flooring #
	###################################
	print "\nTesting nLstrip"
	
	if nLstrip( 'NNNNNNTATGAAGCCTGTGACGGTTAGCTGTCGGAGGAAGTTTACCCTAG', 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWOvJJ]tiJ', 'B') == 'BBBBBBKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWOvJJ]tiJ':
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s nLstrip 1/1" % (currStatus)
	
	###############################
	# Remove 3' low-quality bases #
	###############################
	print "\nTesting trimLowQualBases"
	
	tmpSeq, tmpQual = trimLowQualBases( 'ACGCGATATGAAGCCTGTGACGGTTAGCTGTCGGAGGAAGTTTACCCTAG', 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWO@ABCDEF', qcutoff=70 )
	
	if tmpSeq == 'ACGCGATATGAAGCCTGTGACGGTTAGCTGTCGGAGGAAGTTT' and tmpQual == 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWO':
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s trimLowQualBases 1/1" % (currStatus)
	
	print "%s trimLowQualBases 1/1" % (currStatus)
	
	##############################
	# Trim adapters off sequence #
	##############################
	print "\nTesting eatAdapters"
	
	# adapter 1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
	# adapter 2 AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
	# fxn returns seq, qual, true/false depending whether seqs are 'good' or not
	
	localAdapterList = ('AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA', 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG')
	localReList = [ SearchAdapter(localAdapterList[0], 0, 5), SearchAdapter(localAdapterList[0], 8, 5), SearchAdapter(localAdapterList[1], 0, 5), SearchAdapter(localAdapterList[1], 8, 5) ]
	
	#reList=localReList
	#localMinAdaptMatchLength=5
	#subMatrix=SUBSTITUTION_MATRIX
	#adaptScoreCutoff=5
	#adapterSeqList = localAdapterList
	#minSeqSize = 10
	#minAdptIdentity=0.87 #SET
	
	# def eatAdapters(seq, qual, reList, minAdaptMatchLength=RE_LEN, subMatrix=SUBSTITUTION_MATRIX, adaptScoreCutoff=ADAPT_SCORE_CUTOFF, minAdptIdentity=0.87 )
	
	# no match
	seqTmp, qualTmp, logicalWasTrimmed = eatAdapters('ACGCGATATGAAGCCTGTGACGGTTAGCTGTCGGAGGAAGTTTACCCTAG', 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWOvJJ]tiJ', reList=localReList, minAdaptMatchLength=5, subMatrix=SUBSTITUTION_MATRIX, adaptScoreCutoff=5, minAdptIdentity=0.90)
	
	if seqTmp == 'ACGCGATATGAAGCCTGTGACGGTTAGCTGTCGGAGGAAGTTTACCCTAG' and qualTmp == 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWOvJJ]tiJ' and logicalWasTrimmed == 0:
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s eatAdapters (no match) 1/5" % (currStatus)
	
	# match 1 -- 16 bp
	seqTmp, qualTmp, logicalWasTrimmed = eatAdapters('ACGCGATATGAAGCCTGTGACGGTTAGCTGTCGGAGATCGGAAGAGCGTC', 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWOvJJ]tiJ', reList=localReList, minAdaptMatchLength=5, subMatrix=SUBSTITUTION_MATRIX, adaptScoreCutoff=5, minAdptIdentity=0.90)
	
	if seqTmp == 'ACGCGATATGAAGCCTGTGACGGTTAGCTGTCGG' and qualTmp == 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFm' and logicalWasTrimmed == 1:
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s eatAdapters (match index 1) 2/5" % (currStatus)
	
	# match 2 -- 22 bp
	seqTmp, qualTmp, logicalWasTrimmed = eatAdapters('ACGCGATATGAAGCCTGTGACGGTTAGCAGATCGGAAGAGCGGTTCAGCA', 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWOvJJ]tiJ', reList=localReList, minAdaptMatchLength=5, subMatrix=SUBSTITUTION_MATRIX, adaptScoreCutoff=5, minAdptIdentity=0.90)
	
	if seqTmp == 'ACGCGATATGAAGCCTGTGACGGTTAGC' and qualTmp == 'uuW^YYKYyIMnZVGGD_J}~O~GZuMn' and logicalWasTrimmed == 1:
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
		
	print "%s eatAdapters (match index 2) 3/5" % (currStatus)
	
	# long match -- 33 ndx 2 bp with min seq size 25. Overtrimmed, but fxn does not check trimming
	seqTmp, qualTmp, logicalWasTrimmed = eatAdapters('ACGCGATATGAAGCCTGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG', 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWOvJJ]tiJ', reList=localReList, minAdaptMatchLength=5, subMatrix=SUBSTITUTION_MATRIX, adaptScoreCutoff=5, minAdptIdentity=0.90)
	
	if seqTmp == 'ACGCGATATGAAGCCTG' and qualTmp == 'uuW^YYKYyIMnZVGGD' and logicalWasTrimmed == 1:
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s eatAdapters (long adapter trimming) 4/5" % (currStatus)
	
	# too little identity 85% identity, 20bp, require 90% (3 mismatches)
	seqTmp, qualTmp, logicalWasTrimmed = eatAdapters('ACGCGATATGAAGCCTGTGACGGTTAGCTGAGATCGGAAGACCGGTACTG', 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWOvJJ]tiJ', reList=localReList, minAdaptMatchLength=5, subMatrix=SUBSTITUTION_MATRIX, adaptScoreCutoff=5, minAdptIdentity=0.90)
	
	if seqTmp == 'ACGCGATATGAAGCCTGTGACGGTTAGCTGAGATCGGAAGACCGGTACTG' and qualTmp == 'uuW^YYKYyIMnZVGGD_J}~O~GZuMnL]}WFmlYJMXhZWOvJJ]tiJ' and logicalWasTrimmed == 0:
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s eatAdapters (match, but too little identity 5/5)" % (currStatus)
	
	###################################
	# Overall sequence cleaning suite #
	###################################
	#print "Testing cleanSequence"
	
	#######################################################
	# Fxn returns string True or False depending on input #
	#######################################################
	print "\nTesting boolString"
	
	if boolString(True) == "True":
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s boolString (True) 1/2" % (currStatus)
	
	if boolString(False) == "False":
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s boolString (False) 2/2" % (currStatus)
	
	##############################
	# Is the string a valid int? #
	##############################
	print "\nTesting isIntString"
	
	if isIntString('1'):
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s isIntString(\'1\') 1/2" % (currStatus)
	
	if not isIntString('NotANumber'):
		testsPassed = testsPassed + 1
		currStatus = PASS
	else:
		testsFailed = testsFailed + 1
		currStatus = FAIL
	
	print "%s isIntString (\'NotANumber\') 2/2" % (currStatus)
	
	#############
	# Summarize #
	#############
	totalTests = testsPassed + testsFailed
	print "\nUnit testing summary:"
	print "Passed %s/%s (%.1f%%)" % (testsPassed, totalTests, float(testsPassed)/float( nonZero(totalTests) ) * 100.0)
	print "Failed %s/%s (%.1f%%)" % (testsFailed, totalTests, float(testsFailed)/float( nonZero(totalTests) ) * 100.0)
	
	sys.exit(0)

def getSequenceAmount( totalSeqBp ):
	if totalSeqBp > 1E15:
		ext = "Pbp"
		adjBp = totalSeqBp / 1E15
	elif totalSeqBp > 1E12:
		ext = "Tbp"
		adjBp = totalSeqBp / 1E12
	elif totalSeqBp > 1E9:
		ext= "Gbp"
		adjBp = totalSeqBp / 1E9
	elif totalSeqBp > 1E6:
		ext = "Mbp"
		adjBp = totalSeqBp / 1E6
	elif totalSeqBp > 1E3:
		ext = "Kbp"
		adjBp = totalSeqBp / 1E3
	else:
		ext = "bp"
		adjBp = totalSeqBp
	return SequenceInterpretation( totalSeqBp, adjBp, ext )

#############################
# Run unittest if requested #
#############################
if UNITTEST:
	unittest()
	sys.exit(0)

#####################
# Set quality scale #
#####################
# important: the q cutoff DOES NOT require using the offset with this setup
QUALITYSCALE = setQualScale(INFILE[0])
if QUALITYSCALE == 'illumina1.3':
	TRIMCUTOFF += ILMN_OFFSET
	LOW_SCORE = chr( ILMN1P3_MIN )
elif QUALITYSCALE == 'phred':
	TRIMCUTOFF += SNGR_OFFSET
	LOW_SCORE = chr( PHRED_MIN )
elif QUALITYSCALE == 'illumina1.0':
	errMsg( "Zoinks! We've got old illumina log probability scale here. You're on your own buddy." )

######################
# Set junk extension #
######################
if INTYPE == 'scarf':
	JNKEXT = '.scarf'
elif INTYPE == 'fastq':
	if QUALITYSCALE == 'illumina1.3':
		JNKEXT = ".fq"
	elif QUALITYSCALE == 'phred':
		JNKEXT = ".fqs"

#########################
# print settings to log #
#########################
print "%s v%s" % (SCRIPT_PATH, VERSION)
print "Options set"
print "==========="
print "Input-scale: %s" % (QUALITYSCALE)
print "Ends: %s" % (ENDS)
print "Input files: %s" % (', '.join(INFILE))
print "Output base: %s" % (OUTFILE)
print "Output extension: %s" % (OUTPUTEXTENSION)
print "Output quality scale: %s" % (OUTPUTSCALE)
print "Read 1 - left hardclip: %s" % (R1LTRIM)
if R1LTRIM != 0:
	print "**Warning** 5' hardclip can interfere with mark duplicates. Recommend rerunning with left trimming disabled"
print "Read 1 - right hardclip: %s" % (R1RTRIM)
print "Read 2 - left hardclip: %s" % (R2LTRIM)
if R2LTRIM != 0:
	print "**Warning** 5' hardclip can interfere with mark duplicates. Recommend rerunning with left trimming disabled"
print "Read 2 - right hardclip: %s" % (R2RTRIM)
print "Remove adapter: %s" % (boolString(REM_ADAPT))
print "Adapter sequence(s): %s" % (', '.join(ADAPT_SEQ))
print "Right quality trim: %s" % (boolString(RIGHTQUALTRIM))
sys.stdout.flush() # forces writing immediately. otherwise sits in buffer until the program finishes

###################################################
# Build regular expressions for adapter stripping #
###################################################
for currAdapter in ADAPT_SEQ:
	re_start = 0 - RE_STEP
	for reLoop in range(RE_NUM):
		re_start += RE_STEP
		if re_start + RE_LEN >= len(currAdapter):
			break
		SEARCH_LIST.append( SearchAdapter(currAdapter, re_start, RE_LEN) )

########################
# setup files for junk #
########################
if ENDS == 'paired':
	junk_file1 = OUTFILE + "_1_junk" + JNKEXT
	junk_file2 = OUTFILE + "_2_junk" + JNKEXT
	
	try:
		junkFH1 = open(junk_file1, 'w')
	except:
		errMsg( "Could not open first junk file for PE reads [%s]!" % (junk_file1) )
		
	try:
		junkFH2 = open(junk_file2, 'w')
	except:
		errMsg( "Could not open first junk file for PE reads [%s]!" % (junk_file2) )
	
	junkHandles.append( junkFH1 )
	junkHandles.append( junkFH2 )
		
else:
	errMsg("We're having trouble with the ENDS specified [%s]..." % (ENDS))

##########################
# setup std output files #
##########################
if ENDS == 'paired':
	good_file1 = OUTFILE + "_1_clean" + OUTPUTEXTENSION
	good_file2 = OUTFILE + "_2_clean" + OUTPUTEXTENSION
	
	try:
		goodFH1 = open(good_file1, 'w')
	except:
		errMsg( "Could not open first file for clean reads [%s]!" % (good_file1) )
		
	try:
		goodFH2 = open(good_file2, 'w')
	except:
		errMsg( "Could not open second file for clean reads [%s]!" % (good_file2) )
		
	goodHandles.append( goodFH1 )
	goodHandles.append( goodFH2 )
		
else:
	errMsg("We're having trouble with the ENDS specified [%s]..." % (ENDS))

###################
# start the clock #
###################
analysisStartTime = datetime.now()

##############
# single end #
##############
if ENDS == 'single':
	errMsg( "Single-end mode not currently supported." )

##############
# Paired end #
##############
elif ENDS == 'paired':	
	if INTYPE == 'fastq':
		###############
		# Open input1 #
		###############
		try:
			inputFile1 = open(INFILE[0], 'r')
		except:
			errMsg( "Could not open input file %s" % (INFILE[0]) )
		
		###############
		# Open input2 #
		###############
		try:
			inputFile2 = open(INFILE[1], 'r')
		except:
			errMsg( "Could not open input file %s" % (INFILE[1]) )
		
		######################
		# Parse line by line #
		######################
		while 1:
			
			##########################################################################
			# not being elegant here. big boo boo. Hopefully it won't eat you alive. #
			# don't ever do a while True or while 1.                                 #
			# that said, take your judgement elsewhere.                              #
			##########################################################################
			
			##############
			# read in id #
			##############
			try:
				currentId1 = inputFile1.readline().rstrip('\r\n')
			except:
				break
			
			if len(currentId1) == 0:
				break
			
			####################
			# read in sequence #
			####################
			currentSeq1 = inputFile1.readline().rstrip('\r\n')
			
			if len(currentSeq1) == 0:
				errMsg( "Problem reading sequence line from fastq" )
			
			#####################
			# read in ID line 2 #
			#####################
			repeatId1 = inputFile1.readline().rstrip('\r\n')
			
			if len(repeatId1) == 0:
				errMsg( "Problem reading in second ID line" )
			
			###################
			# read in quality #
			###################
			currentQuality1 = inputFile1.readline().rstrip('\r\n')
			
			if len(currentQuality1) == 0:
				errMsg( "Problem reading in quality" )
			
			##############
			# Now read 2 #
			##############
			##############
			# read in id #
			##############
			currentId2 = inputFile2.readline().rstrip('\r\n')
			
			if len(currentId2) == 0:
				errMsg( "Problem reading read 2 identifier" )
			
			####################
			# read in sequence #
			####################
			currentSeq2 = inputFile2.readline().rstrip('\r\n')
			
			if len(currentSeq2) == 0:
				errMsg( "Problem reading sequence line from fastq" )
			
			#####################
			# read in ID line 2 #
			#####################
			repeatId2 = inputFile2.readline().rstrip('\r\n')
			
			if len(repeatId2) == 0:
				errMsg( "Problem reading in second ID line" )
			
			###################
			# read in quality #
			###################
			currentQuality2 = inputFile2.readline().rstrip('\r\n')
			
			if len(currentQuality2) == 0:
				errMsg( "Problem reading in quality" )
			
			###################
			# Check read sync #
			###################
			read1_id = currentId1.split('#', 1)
			if len(read1_id) == 1:
				read1_id = currentId1.split(' ', 1)[0]
			else:
				read1_id = read1_id[0]
				
			read2_id = currentId2.split('#', 1)
			if len(read2_id) == 1:
				read2_id = currentId2.split(' ', 1)[0]
			else:
				read2_id = read2_id[0]
			
			if read1_id != read2_id:
				errMsg( "Reads out of sync!!!\nRead1: [%s]\nRead2: [%s]\n" % (read1_id, read2_id) )
			
			######################
			# Count read lengths #
			# Iter read number   #
			######################
			readCount += 1
			totalSequenceCount += len(currentSeq1) + len(currentSeq2)
			
			#################################
			# Get N percent from seq length #
			#################################
			n_perc1 = float(currentSeq1.count('N')) / float(len(currentSeq1))
			n_perc2 = float(currentSeq2.count('N')) / float(len(currentSeq2))
			
			if n_perc1 >= N_PERC_CUTOFF and n_perc2 >= N_PERC_CUTOFF:
				# write to junk
				highNsCount += 1
				
				junkHandles[0].write("%s\n%s\n+\n%s\n" % (currentId1, currentSeq1, currentQuality1))
				junkHandles[1].write("%s\n%s\n+\n%s\n" % (currentId2, currentSeq2, currentQuality2))
				
			elif len(currentSeq1) != len(currentQuality1) or len(currentSeq2) != len(currentQuality2):
				escapesOfDeath += 1
				
				junkHandles[0].write("%s\n%s\n+\n%s\n" % (currentId1, currentSeq1, currentQuality1))
				junkHandles[1].write("%s\n%s\n+\n%s\n" % (currentId2, currentSeq2, currentQuality2))
			
			else:
				#################################
				# Clean the data before writing #
				#################################
				cleaned_read1 = cleanSequence( currentSeq1, currentQuality1, "read1", trimVal=TRIMCUTOFF, adaptReList=SEARCH_LIST, lowQualScore=LOW_SCORE )
				adapterTrimmedCount += cleaned_read1.incr_adapt
				
				################
				# Clean read 2 #
				################
				cleaned_read2 = cleanSequence( currentSeq2, currentQuality2, "read2", trimVal=TRIMCUTOFF, adaptReList=SEARCH_LIST, lowQualScore=LOW_SCORE )
				adapterTrimmedCount += cleaned_read2.incr_adapt
				
				#################################################
				# Action depends on if read 1 or read 2 are bad #
				#################################################
				if cleaned_read1.good_seq == False and cleaned_read2.good_seq == False:
					#############################
					# Both too short. Junk 'em. #
					#############################
					badReadCount += 1
					
					junkHandles[0].write("%s\n%s\n+\n%s\n" % (currentId1, currentSeq1, currentQuality1))
					junkHandles[1].write("%s\n%s\n+\n%s\n" % (currentId2, currentSeq2, currentQuality2))
					
					continue
				
				##########################################
				# Convert qualities and junk bad escapes #
				##########################################
				try:
					r1CleanQual = qualityConverter(cleaned_read1.quality, inscale=QUALITYSCALE, outscale=OUTPUTSCALE)
					r2CleanQual = qualityConverter(cleaned_read2.quality, inscale=QUALITYSCALE, outscale=OUTPUTSCALE)
				except KeyError:
					escapesOfDeath += 1
					
					junkHandles[0].write("%s\n%s\n+\n%s\n" % (currentId1, currentSeq1, currentQuality1))
					junkHandles[1].write("%s\n%s\n+\n%s\n" % (currentId2, currentSeq2, currentQuality2))
					
					continue
				
				usableSequenceCount += cleaned_read1.goodBp + cleaned_read2.goodBp
					
				###################################
				# Write line to proper index file #
				###################################
				goodHandles[0].write("%s\n%s\n+\n%s\n" % ( currentId1, cleaned_read1.sequence, cleaned_read1.quality))
				goodHandles[1].write("%s\n%s\n+\n%s\n" % ( currentId2, cleaned_read2.sequence, cleaned_read2.quality))
				
				##########################
				# Confirm index match    #
				# Short-circuit the loop #
				##########################
				match = True
			
		inputFile1.close()	
		inputFile2.close()
		
	else:
		errMsg( "Input type [%s] not recognized" % (INTYPE) )
		
#########################
# Not recognized option #
#########################
else:
	errMsg( "End option [%s] not recognized" % (ENDS) )

#############################
# Print summary information #
#############################
nominalSequence = getSequenceAmount( totalSequenceCount )
usableSequence = getSequenceAmount( usableSequenceCount )

print """
Cleaning completed in %s seconds

###########
# Summary #
###########
""" % (datetime.now() - analysisStartTime)
# Total X end reads processed
print "%s total %s-end reads processed" % (readCount, ENDS)

print "%s reads (%.1f%%) ignored due to high 'N' content" % (highNsCount, float(highNsCount)/float( nonZero(readCount) )*100.0)

print "%s reads (%.1f%%) lost due to escape characters in quality scores" % (escapesOfDeath, float(escapesOfDeath) / float( nonZero(readCount) ) *100.0)

print "%s reads (%.1f%%) dropped due to overtrimming" % (badReadCount, float(badReadCount)/float( nonZero(readCount) )*100.0)

if REM_ADAPT == True:
	print "%s individual reads (%.1f%% of total ends) had adapter sequence stripped" % (adapterTrimmedCount, float(adapterTrimmedCount) / float( nonZero(readCount) ) * 100.0 / (2.0 if ENDS == 'paired' else 1.0))

print "%.1f %s total sequence" % (nominalSequence.adjustBp, nominalSequence.extensionName)

print "%.1f %s usable sequence* (%.1f%%) after cleaning" % (usableSequence.adjustBp, usableSequence.extensionName, float(usableSequenceCount) / float(totalSequenceCount) * 100.0 )

print
print "*All Ns ignored in usable sequence count"

###########
# Cleanup #
###########
for FH in goodHandles:
	FH.close()

for FH in junkHandles:
	FH.close()
