Contents:

The ch.* files are standard variety of data files with various parts
(word stats, PMI, etc.) used in hca.  These files were created using
linkBags and other routines.

The cht.* files are an experimental version to test out the collocation
reporting facilities, which use ".txtbag" files in sequential format
and ".colls" files that report on the candidate collocations in each bag.
These were created with an experimental feature of linkTables and linkBags
  linkTables --mindf 4 --colllastnotstop --collsize 3 --stopfile ~/Desktop/stops.txt church.links ch
  linkBags --nobag --repcoll  church.links ch
  #   121 = #stopwords reported in ch.words
  scripts/bagstop.pl ch cht 121
