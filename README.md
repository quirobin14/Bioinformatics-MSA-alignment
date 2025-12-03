# Bioinformatics-MSA-alignment

This is a comparison between a normal MSA alignment making a whole DP table and an alignment that is looked through using the A* algorithm

(paper about it found at https://cdn.aaai.org/AAAI/2002/AAAI02-111.pdf)

There is a main method printing out a lot of information such as 

-time
-memory usage
-alignments themselves

to compare them to see where each algorithm shines

also to run the code download it to your device and compile it
(javac msaAlignments.java)

then run the command -> 
java msaAlignments.java <input.fasta>

or if you want a custom output file -> 
java msaAlignments <test.fasta> <output.txt>

