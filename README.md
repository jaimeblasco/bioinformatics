# bioinformatics

Bioinformatic Methods I - https://www.coursera.org/learn/bioinformatics-methods-1/

## BLAST
https://blast.ncbi.nlm.nih.gov/Blast.cgi

Program Query Database Alignment # Searches Uses
blastn DNA DNA DNA 1 find homologous DNA sequences
tblastx DNA DNA protein 36 find homologous proteins from unannotated
query and db sequences
blastx DNA protein protein 6 identify proteins in query DNA sequence
tblastn protein DNA protein 6 find homologous proteins in unannotated
DNA DB
blastp protein protein protein 1 find homologous proteins

## MULTIPLE SEQUENCE ALIGNMENT 

Mega X
https://www.megasoftware.net/

Global alignments are those in which the set of sequences to be aligned are similar across their
entire length
Local alignments are those in which there is similarity only at particular sub-regions of the
sequence. 

Clustal
DIALIGN - http://dialign.gobics.de/
MAFFT - http://mafft.cbrc.jp/alignment/server/

## Phylogenetics

MEGA X http://www.megasoftware.net/
PHYLIP http://evolution.genetics.washington.edu/phylip.html
Web-based ML analysis http://bar.utoronto.ca/webphylip/ 

## Selection Analysis

dN/dS Ratio Test
It is often useful to check your aligned DNA sequences using the TRANSEQ app of EMBOSS suite of tools from the EBI. Go to http://bar.utoronto.ca/EMBOSS
Datamonkey (http://www.datamonkey.org/) to look for selection in our sequences

## 'Next Gen' Sequence Analysis (RNA-Seq) / Metagenomics

RNA-Seq Analysis

Plant MPSS databases homepage, at http://mpss.danforthcenter.org/ 1

METAGENOMICS

http://metagenomics.anl.gov/

Bioinformatic Methods II - https://www.coursera.org/learn/bioinformatics-methods-2/

##  Protein Domain, Motif and Profile Analysis 

CDD – Conserved Domain Database 
http://www.ncbi.nlm.nih.gov/

SMART
http://smart.embl-heidelberg.de/

Pfam
http://pfam.xfam.org/

InterproScan
http://www.ebi.ac.uk/interpro/search/sequence-search

## PPI
List PPI databases:
http://ppi.fli-leibniz.de/jcb_ppi_databases.html

DIP 
http://dip.doe-mbi.ucla.edu/dip/Main.cgi

BIOGRID
http://www.thebiogrid.org/index.php

Cytoscape
http://www.cytoscape.org/download.html


Bioinformatic Methods I - https://www.coursera.org/learn/bioinformatics-methods-2/

## Protein Domain, Motif and Profile Analysis

know why we are interested in searching for motifs and profiles in sequences;
 know the advantages and disadvantages of representing structural elements in protein
sequences as motifs or profiles;
 be able to generate a motif given an alignment;
 understand how to score a given sequence with a PSSM or profile HMM;
 be able to use CDD, CDART, SMART, Pfam, and InterProScan to identify specific
functional units within protein sequences; 

CDD – Conserved Domain Database 
http://www.ncbi.nlm.nih.gov/

CDART

SMART
http://smart.embl-heidelberg.de/

Pfam
http://pfam.xfam.org/

InterProScan
http://www.ebi.ac.uk/interpro/search/sequence-search

## Protein-Protein Interactions (PPIs)

understand why protein-protein interactions are important biologically, and also how they
may be determined experimentally;
 be able to assess the advantages and disadvantages of the methods for determining
protein-protein interactions;
 know the terminology associated with protein-protein interaction graphs;
 be able to use DIP, BioGRID and Cytoscape to identify interacting proteins for your gene
product of interest and to filter and decorate networks based on additional information;
 be able to identify the type of support for a given interaction in a given database;
 be able to interpret the other types of information (GO categories) provided by the
software tools. 

DIP
http://dip.doe-mbi.ucla.edu/dip/Main.cgi

BioGRID
http://www.thebiogrid.org/index.php

Cytoscape

## Protein Structure

know the main methods for determining protein structure;
 be familiar with Protein Database records and how to determine which method was used
to ascertain a given protein’s structure;
 be able to view the protein’s structure both with the Jsmol applet and PyMOL;
 be able to use PyMOL to view the structure model in different representations (ribbon,
stick etc.) and colour different parts of the molecule at will;
 be able to align two structures using the command line interface of PyMOL;
 be able to highlight certain residues in the Sequence Viewer part of PyMOL and have
these displayed in the Structure Viewer window
 be able to measure inter-atom distances using the Measure tool of PyMOL. 

Protein DataBank
https://www.rcsb.org/

NCBI Vast
http://www.ncbi.nlm.nih.gov/Structure/VAST/vastsearch.html

PyMOL

##  Gene Expression Analysis 

know the main technologies for creating expression data and for extracting information
from microarray hybridizations or RNA-seq reads;
 understand the sources of error associated with the outputs of a transcriptomics
experiment;
 understand the importance of normalization and how to interpret MA and box plots;
 be able to select significant genes and to organize them using hierarchical clustering –
what does the Pearson correlation coefficient score measure?;
 know how to use gene expression databases to leverage preexisting expression data;
 be able to use BioConductor for normalization, selecting significantly differentially
expressed genes, creating heatmaps, and clustering;
 be able to use some online tools to identify coexpressed genes using existing publiclyavailable gene expression data and to view their expression patterns;
 be familiar with gene expression databases. 

Gene Expression Omnibus (GEO) and the Sequence Read Archive (SRA)
http://www.ncbi.nlm.nih.gov/sra

BioConductor
http://bar.utoronto.ca/BioC/

Expression Angler 
http://bar.utoronto.ca/ntools/cgi-bin/ntools_expression_angler.cgi

Expression Browser
http://bar.utoronto.ca/affydb/cgi-bin/affy_db_exprss_browser_in.cgi

AgriGO
http://systemsbiology.cau.edu.cn/agriGOv2/classification_analysis.php?category=Plant&&family=Brassicaceae

eFP Browser 
http://bar.utoronto.ca/efp/

ATTED-II 
http://atted.jp/

GeneMANIA
http://genemania.org/

## Cis Regulatory Element Mapping and Prediction

know the main technologies for identifying transcription factor binding sites in vitro;
 be able to name some methods for identifying potential transcription factor binding sites
in silico;
 understand how word-count and Gibbs sampling methods work and the differences
between these two methods;
 understand why we would want to generate results which referenced a background set of
all promoters;
 know how to use online tools, such as TF2Network or JASPAR to identify known ciselements in promoters of coexpressed genes;
 be able to use online versions of MEME. 

Promomer
http://bar.utoronto.ca/

TF2Network
http://bioinformatics.psb.ugent.be/webtools/TF2Network/

JASPAR and Cistome 
http://bar.utoronto.ca/cistome_legacy/cgi-bin/BAR_Cistome.cgi
http://jaspar.genereg.net/

MEME
http://meme-suite.org/tools/meme

Human eFP Browser
http://bar.utoronto.ca/efp_human/cgi-bin/efpWeb.cgi


Plant Bioinformatics - https://www.coursera.org/learn/bioinformatics-methods-1/





