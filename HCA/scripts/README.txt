Scripts
-------

All these scripts are offered up as useful, but they are not well supported,
not well documented, and not well integrated.  You will also need to 
modify the header for your particular Perl installation.

The main scripts come with their own man page so use the 
option "--man" to get the manpage:
e.g.,   ./cooc2pmi.pl --man
        perldoc ./cooc2pmi.pl 
Though for this to work, you must have Perl's "pod2usage" 
and "perldoc" programmes installed.

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

Routine used to strip first K words from a vocab.  Modifies following files:
  .words, .tokens, .txtbag, .colls
Note only works with sequential version of .txtbag.
Be warned ...

Routine to generate word clouds and graph.  Requires other code and
only tested in Linux but should work on Mac OSX.
1.  Make sure you have graphviz installed.
2.  Install Andreas Mueller's wordcloud in Python, for instance
    using pip(3).
3.  Replace his "wordcloud.py" library file with the one in this directory.
    Just hope the versions match!
4.  Copy "topset2word.pl" from this directory into your path.
4.  Then your are good.  Read the man page.  See example outputs at:
       http://topicmodels.org/2016/03/25/visualising-a-topic-model/

Others
--------
Some handy routines with a wee bit of documentation
but no self-contained man entries.
    perp.pl
    mkselfbetapr.pl
    pretty.xsl
    drep.pl
    msx2ms.pl
    dtm2tca.pl
    prldac.pl  - does precision-recall calcs for test data based on predictions fromm hca
    trim.pl    - converts output of "hca -v -v -V" topic report to handy form
