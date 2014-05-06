Scripts
-------

The main scripts come with their own man page so use the 
option "--man" to get the manpage:
e.g.,   ./cooc2pmi.pl --man

Generate PMI, NPMI or robust version of PMI from the output of
of linkCoco (see DCA-Bags):
      cooc2pmi.pl

Take the above produced PMI/NPMI file and specialise it
for the current data set (i.e.,  project onto the specific
dictionary):
      mkmat.pl

Process the ".word" output from linkTables to create a 
proportion vector matching the dictionary file.
Can be used with the "-u" option to "hca".
      mkbetapr.pl

Run different tests using hca with
      runtest.pl

When a "-lsp,2,1" style command is used with "hca" it
reads the STEM.smap file for a list of word indices and then
creates sparsity data (i.e. topic strengths) for each of
the words.  This is output to RESTEM.smap.  Create a
nice PNG figure with these using:
       plotsmap.pl

Routines process the sparse matrix format of 
      http://archive.ics.uci.edu/ml/datasets/Bag+of+Words
   spcat.pl -- merge two files
   spcut.pl --  take first N records out of file
   spprint.pl -- print a matrix with words, assuming both rows and columns
                are indexes to the vocab file
   spthin.pl -- thin out vocabulary

Routine used with collections having an ".epoch", to generate
a train-test split and create the matching ".epoch" file.
Works with ldac or txtbag data only.
   berand.pl


Others
--------
Some handy routines with a wee bit of documentation
but no self-contained man entries.
    perp.pl
    mkselfbetapr.pl
    pretty.xsl
    drep.pl
    msx2ms.pl
