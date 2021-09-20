### Processing FASTA files

Program processing FASTA files can be found in toolbox/fasta.py file. It is in a form of commandline app accepting following parameters:

`-f / --file` for specifying the input file\
`-d / --desc` for printing the description of given molecule\
`-q / --seq` for printing the sequence of given molecule\
`-l / --length` for printing the length of the sequence of given molecule\
`-s / --subseq` for printing a subsequence of given sequence defined by start and end index\


the `-f` option is obligatory for each call, the rest is optional\
After calling, it automatically parses the FASTA file and runs appropriate functions.

Examples:
```commandline
fasta.py -f "data/sequence" -q "QGKT01000001.1"
fasta.py -f "data/sequence" -s "QGKT01000001.1" 34 789
```


### Measuring sequence similarity using Hamming distance

Program measuring sequence similarity can be found in toolbox/hamming.py file. It is also a commandline app. It accepts the following parameters.\
<br>
`-f / --file` for specifying the input file\
`-s / --seqs` for specifying the measured sequences

Example:
```commandline
hamming.py -f "data/sequence" -s "TEST1" "TEST2"
```
I have modified the "sequence" file by adding two made up sequences named "TEST1" and "TEST2" of the same length for the testing sake.


### Sequence alignment using edit distance

Program doing sequence alignment can be found in toolbox/SA.py file. This file is no longer a command app but a library.\
The sequence alignment is done by function align which takes two sequences as strings ands prints all sequence alignments with the same minimum score.
```python
import SA
SA.align("CLOCK", "LACKS")
```
### Processing PDB files

Library for processing PDB files can be found in toolbox/pdb.py file. 
The main object in this library is Structure. It represents the whole input PDB file.
```python
import pdb
strctr = pdb.Structure("A2a.pdb")
```
The PDB file is represented in two different ways. One is with a respect to structure that roots in the Structure object and branches respectively into models, chains, residues and atoms. The other works like a summary of all the models, chains, residues and atoms included in the PDB file. Those are all listed directly in the Structure object as its attributes. Attribute models is included in both of those approaches. It contains all models from the file and also serves as nodes of the tree structure. 
```python
import pdb
strctr = pdb.Structure("A2a.pdb")
models = strctr.models
chains = strctr.chains
residues = strctr.residues
atoms = strctr.residues
#each variable contains a list of instances of corresponding classes. 
#The instance made by PDBParser may be obtained as follows
atom = strctr.atoms[0].atom
```
Then we can compute the width of the structure
```python
import pdb
strctr = pdb.Structure("A2a.pdb")
width = strctr.width
```
And obtain list of atoms or residues being in given distance from given ligand. The desired element is given to the function in the last argument "wanted" as "A" for atoms or "R" for residues.
```python
import pdb
strctr = pdb.Structure("A2a.pdb")
nearAtoms = strctr.search_neighbors(strctr.atoms[6], radius=5, wanted="A")
```
### Processing multiple sequence alignment
Processing multiple sequence alignment can be done with the help of the library in toolbox/SA.py.We can read MSA file by giving it to the constructor of Alignment object.
We can also retrieve sequence by its position or ID, retrieve given column from the MSA or retrieve sum of pairs score of a column and whole MSA with respect to given scoring matrix.
the last two functions require stating the type of scoring matrix according to the Bio.Align.substitution_matrices package and penalty used for evaluation of insertion or deletion.
```python
import SA
alignment = SA.Alignment("data/p53_clustal")
seqByID = alignment.getSeqID("UniRef90_C0PUM1")
seqByPosition = alignment.getSeqPos(6)
sumOfPairsCol = alignment.getPairScoreCol(15, "BLOSUM62", 0)
sumOfPairsMSA = alignment.getMSAScore("BLOSUM62", 0)
```
### Conservation determination from multiple aligned sequences
This problem can be solved using the toolbox/SA.py.
we can use it for obtaining a conservation score for a single position, or getting top n most scoring position in the whole MSA. the conservation score is computed by dividing the number of occurrences of the most abundant amino acid by the number of aligned sequences
```python
import SA
alignment = SA.Alignment("data/p53_clustal")
SPscore = alignment.getConservationScoreOnPos(20)
topN = alignment.getTopNConserved(10)
```

### Computing structure-related properties
Computing structure-related properties can be done by the functions in the toolbox/pdb.py
We can compute the ratio of surface and buried amino acids, 
the ratio of amino acids in the core and on the surface of the protein,
and amino acids' composition of buried and exposed amino acids.
```python
import pdb
strctr = pdb.Structure("A2a.pdb")
surfBurDict = strctr.getSurfaceBuriedRatio()
ratio = surfBurDict["buried"]["ratio"]
buriedComposition = surfBurDict["buried"]["distribution"]
surfaceComposition = surfBurDict["surface"]["distribution"]
surfacePolar = surfBurDict["surface"]["surface polar"]
buriedPolar = surfBurDict["buried"]["buried polar"]
```
#### A2a vs hemoglobin
Using the previous function we got following data.\
In the A2a receptor 56% of all amino acids are buried, 29% of buried amino acids are polar and 49% of surface amino acids are polar.
In the hemoglobin protein 45% of all amino acids are buried, 28% of buried amino acids are polar and 49% of surface amino acids are polar.
Proportions of polar to nonpolar amino acids on the surface and inside are in both proteins roughly the same, but in the ratio of all buried aa to surface aa
they do slightly differ. Hemoglobin has more buried aa than A2a receptor, therefore is very likely more compact