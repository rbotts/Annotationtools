from Bio.Seq import Seq, MutableSeq, translate
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
#from Bio.SeqUtils import GC
from Bio import SeqIO, Entrez, Alphabet
#from Bio.Graphics import GenomeDiagram
#from Bio.Graphics.GenomeDiagram import CrossLink
#from reportlab.lib.colors import *
#from reportlab.lib import colors
#import urllib2
import os
import csv
#import numpy as np
#import datetime

#from user_plasmids import *
#from expanding import *
#from copy import deepcopy

###########################################


# routines for sequence annotation


def add_features_fn(seqname):
	# open fasta file with a single contig of interest (.fasta) in fasta format and
    # associated file of nucleotide sequences of all features saved in multi-fasta format (.fnn)
    # (output from GenemarkS)
	# 
	# inserts contents of fnn file as seq features on the seq object
	# generates an output genbank formatted file of the same name with features to be used later
	

    sqname = seqname
	
	# read in seq to add features to
    handle = open(sqname+".fasta","rU")
    record = SeqIO.read(handle,"fasta")
    handle.close()

    record.id = sqname

	#read contents of gmm file.  gmm output is not standard tab or comma delimited file
    gmm_handle = open(sqname+".fnn", "rU")

    features_rec = list(SeqIO.parse(gmm_handle, "fasta"))
    for fture in features_rec:
        feat_info = fture.id.split('|')
        if feat_info[-3] == '+':
            feat_dir = 1
            start = int(feat_info[-2])-1
            stop = int(feat_info[-1])
					
        else:
            feat_dir = -1
            start = int(feat_info[-2])-1
            stop = int(feat_info[-1])
					
        record.features.append(SeqFeature(id=feat_info[0],location=FeatureLocation(start,\
							stop,strand=feat_dir)))
	
    print "Extracted %i features" % len(record.features)
    return record

def Blast_Features(seqname):
	# open contents of seq file, add features from file of protein orfs produced by GeneMarkS
	# blast each feature and returns a single output table of best BLAST hits
	#
    # Currently loads GeneMarkS ORFs onto a record and then blasts those from the record.
    # This step was added so that it was possible to BLAST the translated or untranslated version
    # For efficiency, this step could be removed.
    #
	# Inputs:
    #      seqname.fasta - fasta file holding the entire sequence for annotation
    #      seqname.fa  -
    #
    # Outputs:
    #      seqname.tab - a tab delimited table containing the best BLAST alignments for each ORF
    
    import datetime
    max_hits = 5
    blast_type="blastp" # set type of blast: n, p or x, translate seq record for blastp
    record=add_features_fn(seqname)
    outFileName=seqname+'.tab'
    geneList=[]
    outfile = open(outFileName, "a") # use a to avoid overwriting old data
    outfile.write(seqname+'\t'+datetime.datetime.now().strftime("%m/%d/%y")+'\n')
    outfile.write('ORF\t'+'start\t'+'stop\t'+'direction\t'+'size (bp)\t'+'size (aa)\t'+'Gene\t'+\
		'definition\t'+'GB \t'+'Type \t'+'query cov\t'+'query id\t'+'comment\n')

    for feature in record.features:
		#call blast with nr database
        # print feature.id
        print 'BLAST-ing ' + feature.id
        print feature.extract(record).seq.translate(table="Bacterial", cds=True)
		# USe the first for blast n
		#result_handle = NCBIWWW.qblast(blast_type, "nr", feature.extract(record).seq)
		# use for blastp
        # the taxid limits the search to Enterobacteriaceae group (taxid:91347)
        result_handle = NCBIWWW.qblast(blast_type, "nr", feature.extract(record).seq.translate(table="Bacterial", cds=True), alignments=max_hits)#, entrez_query = '[taxid]91347')

        giList = []
    
        #open file for writing and then save the blast results
        #xml format is chosen because parser is more stable
        save_file = open (seqname+"temp_blast.xml", "w")
        save_file.write(result_handle.read())
        save_file.close()
        result_handle.close()
		
        # open contents of Blast output and write into nice tabular format
        # works the same as Brook_Blast, but output format is different
        # outfield = open(outFileName, "w") #open for appending so don't write over previous results

        #bring file contents back into handle for reading
        result_handle = open(seqname+"temp_blast.xml")
        blast_record = NCBIXML.read(result_handle)
        E_VALUE_THRESH = 0.04
        found_best = False
        count = 0

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    #Extract accession number and description from title
                    title_info = alignment.title.split('|')
                    accession=title_info[3]
                    descript=title_info[4]
                    gi_num = title_info[1]
                    print

                    #find the best match that is not hypothetical
                    #save all others with same max_iden and query cover as that best
                    if not (gi_num in giList): # and descript.find('unknown')==-1 and descript.find('hypothetical')==-1:
                        giList.append(gi_num)
                        #calc max iden and query cover
                        cover = (len(hsp.query)*1.0/len(feature))
                        max_iden = (hsp.identities * 1.0)/len(hsp.query)

                        #save the information about the best match
                        if not found_best:
                            best_cover=cover
                            best_iden=max_iden
                            found_best=True
                            print ("query cover: %s    max_iden: %s" % (cover, max_iden))

                        #collect all other alignments that have same "best" scores
                        #if cover==best_cover and max_iden==best_iden:
                        # keep the top 10 best alignments
                        if len(giList) < max_hits+1:
                                #calc subject start and end based on strand directions
                            if hsp.sbjct_start < hsp.sbjct_end:
                                sbjct_end = hsp.sbjct_end;
                                sbjct_start = hsp.sbjct_start
                            else:
                                sbjct_end = hsp.sbjct_start
                                sbjct_start = hsp.sbjct_end
                            geneList.append((gi_num,sbjct_start,sbjct_end, accession, descript, cover, max_iden))
                            print 'Hit to '+ descript
                            count=count+1
                            outfile.write("%s\t %i\t %i\t %i\t %i\t %i\t \t %s\t %s\t \t %f\t %f\t \n" \
                            % (feature.id, int(feature.location.start), int(feature.location.end),\
                            feature.strand, len(feature), (len(feature))/3, descript, accession,\
                            cover, max_iden))
                            print("%s \t %s \t %s \t %s \t %s" % (gi_num,accession,descript,\
								str(sbjct_start),str(sbjct_end)))
        # if no suitable alignments are  found insert a blank line in the annotation  table
        if giList==[]:
            outfile.write("\n")
        #else:
        #    outfile.write("%s\t %i\t %i\t %i\t %i\t %i\t \t \t \t \t \t \t \n" \
        #                     % (feature.id, int(feature.location.start), int(feature.location.end),\
        #                     feature.strand, len(feature), len(feature)/3))
			
         #   outfile.write("\n") #print a space between orfs
        print("----------------------- ORF BLAST COMPLETE ------- %s records found----------" % (str(count)))
        #outfile.write("\n\n")
    outfile.close()

def Blast_Features_Local(seqname):
	# open contents of seq file, add features from file of orfs
	# blast each file and return a single output table
	#
	# loop through ORFs and print results in outFile
    import datetime

    outcsv=True #True for outputting to csv file, false produces tab delimited output
    blast_type="blastp" # set type of blast: n, p or x, translate seq record for blastp
    blast_db = "nr" # use nr for protein alignments and blastx, use nt for blastn
    blast_db_path = "/usr/local/ncbi-blast-2.2.29+/db/"
    min_evalue = 0.0001
    max_hits = 5
    
    print 'Running %s, using min e-value of %i, and retrieving a maximum of %i hits for each feature' % (blast_type, min_evalue, max_hits)
    
    record=add_features_fn(seqname)
    outFileName=seqname+'.csv'
    geneList=[]
    outfile = open(outFileName, "a") # use a to avoid overwriting old data
    outfile.write(seqname+'\t'+datetime.datetime.now().strftime("%m/%d/%y")+'\n')
    outfile.write('ORF\t'+'start\t'+'stop\t'+'direction\t'+'size (bp)\t'+'size (aa)\t'+'Gene\t'+\
		'definition\t'+'GB \t'+'Type \t'+'query cov\t'+'query id\t'+'comment\n')

    for feature in record.features:
		# call blast with nr database
        print 'BLAST-ing ' + feature.id
        sq = feature.extract(record)
        sq.seq = sq.seq.translate(table="Bacterial", cds=True)
        #sq = SeqRecord(seq=sq, id =feature.id, name =feature.id)
        sq.id=feature.id
        print sq
        #create a temporary fasta file holding only the record for query
		
        SeqIO.write(sq, open('Temp.fasta','w'),'fasta')
        
        #make a call to blast
        # create command line call to blast
        blast_cline = NcbiblastpCommandline(query='Temp.fasta', db=blast_db_path+blast_db, \
								evalue=min_evalue, outfmt=5, out="temp_blast.xml",\
											 max_target_seqs = max_hits)
        
        os.system(str(blast_cline))

        giList = []

		
        # open contents of Blast output and write into nice tabular format
        # output should be manually curated for best alignments
        # outfield = open(outFileName, "w") #open for appending so don't write over previous results

        #bring file contents back into handle for reading
        result_handle = open("temp_blast.xml")
        blast_record = NCBIXML.read(result_handle)
        E_VALUE_THRESH = 0.01
        found_best = False
        count = 0

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    #Extract accession number and description from title
                    title_info = alignment.title.split('|')
                    accession=title_info[3]
                    descript=title_info[4]
                    gi_num = title_info[1]

                    #find the best match that is not hypothetical
                    #save all others with same max_iden and query cover as that best
                    if descript.find('unknown')==-1 and descript.find('hypothetical')==-1 and not (gi_num in giList):
                        giList.append(gi_num)
                        #calc max iden and query cover
                        cover = (len(hsp.query)*1.0/len(feature))
                        max_iden = (hsp.identities * 1.0)/len(hsp.query)

                        #save the information about the best match
                        if not found_best:
                            best_cover=cover
                            best_iden=max_iden
                            found_best=True
                            print ("query cover: %s    max_iden: %s" % (cover, max_iden))

                        #collect all other alignments that have same "best" scores
                        if cover==best_cover and max_iden==best_iden:
                            #calc subject start and end based on strand directions
                            if hsp.sbjct_start < hsp.sbjct_end:
                                sbjct_end = hsp.sbjct_end;
                                sbjct_start = hsp.sbjct_start
                            else:
                                sbjct_end = hsp.sbjct_start
                                sbjct_start = hsp.sbjct_end
                            geneList.append((gi_num,sbjct_start,sbjct_end, accession, descript, cover, max_iden))
                            count=count+1
                            outfile.write("%s\t %i\t %i\t %i\t %i\t %i\t \t %s\t %s\t \t %f\t %f\t \n" \
                            % (feature.id, int(feature.location.start), int(feature.location.end),\
                            feature.strand, len(feature)+1, (len(feature)+1)/3, descript, accession,\
                            cover, max_iden))
                            print("%s \t %s \t %s \t %s \t %s" % (gi_num,accession,descript,\
								str(sbjct_start),str(sbjct_end)))
        # if no suitable alignments are  found insert a blank line in the annotation  table
        if giList==[]:
            outfile.write("\n")
        else:
            if outcsv:
                outfile.write("%s\t %i\t %i\t %i\t %i\t %i\t \t \t \t \t \t \t \n" \
                             % (feature.id, int(feature.location.start), int(feature.location.end),\
                             feature.strand, len(feature)+1, (len(feature)+1)/3))
            else:
                outfile.write("%s\t %i\t %i\t %i\t %i\t %i\t \t \t \t \t \t \t \n" \
						% (feature.id, int(feature.location.start), int(feature.location.end),\
						feature.strand, len(feature)+1, (len(feature)+1)/3))
			
            outfile.write("\n") #print a space between orfs
    print("----------------------- ORF BLAST COMPLETE ------- %s records found----------" % (str(count)))
    outfile.write("\n\n")
    outfile.close()

def annotate_from_BLAST_csv(seqname):
    # Uses csv BLAST output from multiple sequence search and moves hits into a file formatted in the
    # format as BLAST_features method from above
    # assumes BLAST output is named seqname_BLAST.csv
    # In order to run BLAST on a large number of sequences using the protein searh, set taxa Enterobactericiae
    Entrez.email = "rbotts@pointloma.edu"
    
    handle = open(seqname+'BLAST.csv','rU')
    out_handle = open(seqname+'.csv','a')
    outwriter = csv.writer(out_handle, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
    outwriter.writerow([seqname,datetime.datetime.now().strftime("%m/%d/%y")])
    outwriter.writerow(['ORF','start','stop','direction','size (bp)', 'size (aa)', 'Gene',\
        'definition','gb acc', 'Function','query cov', 'query id', 'comment'])
    #out_handle.write(seqname+'\t'+datetime.datetime.now().strftime("%m/%d/%y")+'\n')
    #out_handle.write('ORF\t'+'start\t'+'stop\t'+'direction\t'+'size (bp)\t'+'size (aa)\t'+'Gene\t'+\
	#	'definition\t'+'GB \t'+'Type \t'+'query cov\t'+'query id\t'+'comment\n')
    # read in data
    
    csv_file = csv.reader(handle,delimiter=',',dialect=csv.excel)
    
    for row in csv_file:
    # get seq name, location, length and amino acid length from first entry
        try:
            values = row[0].split('|')
            if values[3] == '+':
                strand = 1
            elif values[3] == '-':
                strand = -1
            accession = row[1].split('|')[3]
            print "Fetching "+accession
            seqhandle=Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
            seq = SeqIO.read(seqhandle,"genbank")
            descript = seq.description

            ln_orf = int(values[5]) - int(values[4]) + 1
            outwriter.writerow([values[0], int(values[4])-1, int(values[5]),\
                strand, ln_orf, ln_orf/3,'',descript, accession, '',\
                float(row[4])*100./(ln_orf/3), float(row[2])])
            #out_handle.write("%s\t %i\t %i\t %i\t %i\t %i\t \t %s\t %s\t \t %f\t %f\t \n" \
                #% (values[0], int(values[4])-1, int(values[5]),\
                #strand, ln_orf, ln_orf/3, descript, accession,\
                #100.*float(row[4])/(ln_orf/3), float(row[2])))
        except (ValueError, IndexError):
            print 'row in table skipped'
    out_handle.close()
    handle.close()


def extract_features_csv(seqname,writegb=True):
	# Extract the annotation features from a csv file of the format output from BLAST features
	#
	# writegb switches writing a genbank record with all of the annotations as an outpur
	# returns a seqrecord with the annotations stored in seqname.csv
    #
    #
    ########  Note that this tool does not work directly on the tab formatted output from BLAST_Features
	


	# open seq record file
	handle = open(seqname+".fasta","rU") 
	record = SeqIO.read(handle,"fasta") 
	handle.close()
	
	record.seq.alphabet=Alphabet.generic_dna #set alphabet for use later
	
	# open annotation file, first two lines are headers 
	ann_handle = open(seqname+".csv","rU")

	# read in data
	csv_file = csv.reader(ann_handle,delimiter=',',dialect=csv.excel)
	
	# extract info in each row and add create temp seq objects to write to image
	record.name=seqname[0:5] # seq record prints the name as locus_id, and there is a maximum length
	record.id = csv_file.next()[0] # skip header row
	print record.id
	csv_file.next() # skip date row
	orfcount=1  # count unnamed orfs for naming
	for row in csv_file:
        
		# catch cases where row contains invalid info 
		print row
		try:
			if row[6]=='':
				feat_name='orf'+str(orfcount)
				orfcount=orfcount+1
			else:
				feat_name = row[6]
			print feat_name
			feat_start = int(row[1])
			feat_end = int(row[2])
			feat_dir = int(row[3])
        	#feat_acc = row[8] Not workingright
			feat_note=row[7]
			if row[9]== '':
				feat_type='CDS'
			else:
				feat_type = row[9]
			record.features.append(SeqFeature(id=feat_name,type=feat_type,location=FeatureLocation(feat_start-1,feat_end,strand = feat_dir),\
                        qualifiers={"gene":feat_name}))

		except ValueError:
			print "A row contains invalid entries, row skipped"

	ann_handle.close()
	for f in record.features:
		print f
	if writegb:
		outhandle=open(seqname+".gb","w")
		SeqIO.write(record,outhandle,"genbank")
		outhandle.close()
		print "Successfully wrote " + seqname + ".gb"
	return record


#   Step by step calls for annotating a sequence.
#   Inputs:
#       name.fasta - fasta sequence of a single sequence to be annotated
#   Workflow:
#       1.  Find likely open reading frames (ORFs) using GeneMarkS.  Save nucleotide (protein sequence) output in fasta format
#           Save file according to convention name.fnn
#       2.  BLAST features.   May use local or web based BLAST (local was too slow)
#                BLAST_features or Blast_Features_Local
#       3.  BLAST tools produces a tab delimited table of the hits.  A user should curate this table to identify which hits make the most sense.
#           When the user is done there should be one hit for each ORF.
#####       Currently the output of step 2 must be changed during step 3 in order to be used for step 4.  This should be fixed
#
#       4.  Upload the final annotation onto the sequence and produce a genbank formatted file with all of the annotations
#           using extract_features_csv

sqname = "TestSequence"

#Blast_Features(sqname)
extract_features_csv(sqname)


