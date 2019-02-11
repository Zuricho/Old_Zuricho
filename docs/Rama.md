# Calculating Ramachandran (phi/psi) Angles

The are at least three python libraries which can be used to load PDB files and calculate the protein backbone's ϕ/ψ angles:

- Konrad Hinsen's [Molecular Modelling Toolkit (MMTK)](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/#MMTK)
  Fussy about loading certain flawed PDB files, but getting phi and psi is very easy.
- Thomas Hamelryck's [Bio.PDB module in BioPython](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/#BioPython)
  Tolerant of odd PDB files, making it a safer bet.
- [Python Macromolecular Library (mmLib)](http://pymmlib.sourceforge.net/)
  Might be worth a look - but I haven't had time.

Provided your PDB files has nothing to funny in it, then I think using MMTK is the easiest way to get protein backbone deihedral angles in python. However, as Bio.PDB is much more tolerant of real world PDB files, it is a better choice for "data mining" tasks like calculating Ramachandran angles.

For some examples of when things go wrong, see the ["Top 500" PDB files](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/top500/). See also this page on [other tools](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/other/).



## Phi/Psi using the Molecular Modelling Toolkit (MMTK)

The following short piece of python code uses Konrad Hinsen's python [Molecular Modelling Toolkit (MMTK)](http://starship.python.net/crew/hinsen/MMTK/index.html) to load a PDB file ([1HMP](http://www.rcsb.org/pdb/explore.do?structureId=1HMP)), calculate ϕ and ψ, and print them to screen (in radians). The code is so succinct because MMTK has support to calculate the psi/phi protein backbone dihedral angles built in:

```
import MMTK.Proteins
protein = MMTK.Proteins.Protein("1HMP.pdb", model="no_hydrogens")
for chain in protein : print chain.name, chain.phiPsi()
```

Or, you can do this residue by residue - printing out a little bit more information as we go:

```
import MMTK.Proteins
protein = MMTK.Proteins.Protein("1HMP.pdb", model="no_hydrogens")
for chain in protein :
    print "%s length %i" % (chain.name, len(chain)),
    print "from %s to %s" % (chain[0].name, chain[-1].name)
    for residue in chain :
        print residue.name, residue.phiPsi()
```

The resulting output will look something like this (depending on which PDB file you use):

```
Warning: Some atoms in a protein have undefined positions.
chain0 length 214 from Ser4 to Ala217
Ser4 (None, 2.3639888877052697)
Pro5 (-1.6218793521220176, 0.22586849933785014)
Gly6 (1.148370920424036, -2.8314425746702248)
...
Lys216 (-2.2379259217992624, -3.10736393340063)
Ala217 (-1.420198039722629, None)
chain1 length 209 from Ser4 to Ala217
Ser4 (None, 2.1932192362631189)
Pro5 (-1.1606220655879966, -0.19493770402923941)
Gly6 (1.4256499097182802, -2.869341157938738)
...
Lys216 (-2.9458191504566242, -2.0288444613318499)
Ala217 (1.3986135690663808, None)
```

Note that for the first residue of each protein chain, phi is undefined, while for the last residue, psi is undefined. Also note that MMTK returns the angles in radians, while most plots are drawn from -180° to +180°.

This simplicity comes at a cost - MMTK is very particular about its PDB files, and will not tolerate much deviation from the standard. In fact, this example 1HMP is a case in point - there is actually a break in Chain B, as we will discover later using [BioPython](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/#BioPython)...

For another example, many PDB files contain "funny" hydrogens (e.g. for undeclared protonated histidines) which will make MMTK choke. One way round this is to explicitly ignore any hydrogens in the PDB file:

```
import MMTK.PDB
import MMTK.Proteins

configuration = MMTK.PDB.PDBConfiguration("1HMP.pdb")
configuration.deleteHydrogens()
protein = MMTK.Proteins.Protein(configuration.createPeptideChains(
                                model = "no_hydrogens"))
...
```

Also, in my experience, MMTK copes far better with files downloaded from the PDB, than it does with those which have been edited by another program. Howver, even files direct from the PDB contain errors and oddities ([more](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/top500/)). As an alternative to MMTK, we could use [BioPython](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/#BioPython).



## Phi/Psi and Residue Type using MMTK

If you know anything about Ramachandran plots, you'll know that the residues are normally sorted into four separate types:

- Glycine (the small side chain makes the protein backbone very flexible
- Proline (their large side chain restricts backbone movement)
- Pre-Proline (proline even messed up any residue before it)
- General (not Proline, not Glycine, not before a Proline)

I used these three functions to categorise the residues:

```
def next_residue(residue) :
    """Expects an MMTK residue, returns the next residue
    in the chain, or None"""
    #Proteins go N terminal --> C terminal
    #The next reside is bonded to the C of this atom...
    for a in residue.peptide.C.bondedTo():
        if a.parent.parent != residue:
            return a.parent.parent
    return None

def residue_amino(residue) :
    """Expects an MMTK residue, returns the three
    letter amino acid code in upper case"""
    if residue :
        return residue.name[0:3].upper()
    else :
        return None

def residue_ramachandran_type(residue) :
    """Expects an MMTK residue, returns ramachandran 'type'
    (General, Glycine, Proline or Pre-Pro)"""
    if residue_amino(residue)=="GLY" :
        return "Glycine"
    elif residue_amino(residue)=="PRO" :
        return "Proline"
    elif residue_amino(next_residue(residue))=="PRO" :
        #exlcudes those that are Pro or Gly
        return "Pre-Pro"
    else :
        return "General"
```

We can then load the PDB file and save the angles as a tab separated file:

```
pdb_code = "1HMP"
output_file = open("%s_mmtk.tsv" % pdb_code,"w")
import MMTK.Proteins
protein = MMTK.Proteins.Protein("%s.pdb" % pdb_code, model="no_hydrogens")
for chain in protein :
    for residue in chain :
        phi, psi = residue.phiPsi()
        if phi and psi :
            #Don't write output when missing an angle
            output_file.write("%s:%s:%s\t%f\t%f\t%s\n" \
                % (pdb_code, chain.name, residue.name,
                   degrees(phi), degrees(psi),
                   residue_ramachandran_type(residue)))
output_file.close()
```

This gives something like this...

```
1HMP:chain0:Pro5	-92.926842	12.941312	Proline
1HMP:chain0:Gly6	65.796807	-162.229709	Glycine
1HMP:chain0:Val7	-81.132082	121.413022	General
...
1HMP:chain1:Lys216	-168.783005	-116.244225	General
```

The `degrees()` function is desribed [below](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/#angles).

If you really want them, you can [download this python script and the sample output](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/#download).

Next: Using these numbers to [draw the Ramachandran Plot with python](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/drawing/) (or [draw it with R](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/ramachandran/)).



## Phi/Psi using BioPython

[BioPython](http://www.biopython.org/) has a PDB parser, which includes its own vector classes which are used internally by the `Polypeptide`class to calculate a protein backbone's torsion angles. My first attempt looked like this:

```
import Bio.PDB
for model in Bio.PDB.PDBParser().get_structure("1HMP", "1HMP.pdb") :
    for chain in model :
        poly = Bio.PDB.Polypeptide.Polypeptide(chain)
        print "Model %s Chain %s" % (str(model.id), str(chain.id)),
        print poly.get_phi_psi_list()
```

The less than perfect output, abbreviated:

```
Model 0 Chain A
[(None, 2.3639878340439484), ..., (-1.4201969110624371, None)]
Model 0 Chain B
[(None, 2.1932190074065052), ..., (1.3986135567166249, None)]
Model 0 Chain  
[(None, None), ..., (None, None)]
```

This "worked" to an extent, but gave a extra polypeptide chain made out of the waters etc. The problem is the that code (above) ignored all the extra validation gained doing it using the `PPBuilder` or `CaPPBuilder` classes:

```
import Bio.PDB
for model in Bio.PDB.PDBParser().get_structure("1HMP", "1HMP.pdb") :
    for chain in model :
        polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
        for poly_index, poly in enumerate(polypeptides) :
            print "Model %s Chain %s" % (str(model.id), str(chain.id)),
            print "(part %i of %i)" % (poly_index+1, len(polypeptides)),
            print "length %i" % (len(poly)),
            print "from %s%i" % (poly[0].resname, poly[0].id[1]),
            print "to %s%i" % (poly[-1].resname, poly[-1].id[1])
            print poly.get_phi_psi_list()
```

Giving the following (which I have abbreviated):

```
Model 0 Chain A (part 1 of 1) length 214 from SER4 to ALA217
[(None, 2.3639878340439484), ..., (-1.4201969110624371, None)]
Model 0 Chain B (part 1 of 2) length 101 from SER4 to TYR104
[(None, 2.1932190074065052), ..., (-1.019038394289536, None)]
Model 0 Chain B (part 2 of 2) length 108 from SER109 to ALA217
[(None, -3.0761010286975696), ..., (1.3986135567166249, None)]
```

Note that this time, the "third chain" (the water molecules etc) does not get turned into any polypeptides.

More interestingly, "Chain B" gets turned into *two* smaller polypeptides... the reason for this is that the `PPBuilder`has checked all the N-C bond lengths, and noticed a big jump between TYR104 and SER109. You'll notice there is a jump in the index number too - consulting the header of the PDB file we see:

A loop of residues 103 - 121 in both chains A and B is poorly ordered. Coordinates given for this region result from a tentative fitting to poor electron density and should be treated with caution. For this loop in the second monomer, residues 105 - 108 and 121 are missing. Some residues in this region are modeled as alanine residues.

So, there really is a jump in Chain B, and TYR104 and SER109 are not bonded together! Something MMTK missed... maybe "Chain B" should have been recorded as two chains to avoid this implied bond?

Anyway, we can do this residue by residue like this:

```
import Bio.PDB
for model in Bio.PDB.PDBParser().get_structure("1HMP", "1HMP.pdb") :
    for chain in model :
        polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
        for poly_index, poly in enumerate(polypeptides) :
            print "Model %s Chain %s" % (str(model.id), str(chain.id)),
            print "(part %i of %i)" % (poly_index+1, len(polypeptides)),
            print "length %i" % (len(poly)),
            print "from %s%i" % (poly[0].resname, poly[0].id[1]),
            print "to %s%i" % (poly[-1].resname, poly[-1].id[1])
            phi_psi = poly.get_phi_psi_list()
            for res_index, residue in enumerate(poly) :
                res_name = "%s%i" % (residue.resname, residue.id[1])
                print res_name, phi_psi[res_index]
```

This gives almost the same results as MMTK. There are minor differences in the floating point values, and the formatting of the identifiers. The only significant difference is BioPython does not give a bogus set of angles for the Tyr104 to Ser109 "bond":

```
Model 0 Chain A (part 1 of 1) length 214 from SER4 to ALA217
SER4 (None, 2.3639878340439484)
PRO5 (-1.621876195380803, 0.225868591662406)
...
LYS216 (-2.2379247741249819, -3.1073635839297786)
ALA217 (-1.4201969110624371, None)
Model 0 Chain B (part 1 of 2) length 101 from SER4 to TYR104
SER4 (None, 2.1932190074065052)
PRO5 (-1.1606194278511266, -0.19494016733756098)
...
SER103 (1.8257513984522575, 0.70759204280141397)
TYR104 (-1.019038394289536, None)
Model 0 Chain B (part 2 of 2) length 108 from SER109 to ALA217
SER109 (None, -3.0761010286975696)
THR110 (1.3981586137324507, -2.0215751441392009)
...
LYS216 (-2.9458182040496657, -2.0288443849169973)
ALA217 (1.3986135567166249, None)
```

We can adapt this code to determine the Ramachandran plot "residue type" and produce a file as done with [MMTK](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/#MMTK)(above). You can [download this code](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/#download) and its output.



## Radians to Degrees

In the sections above we showed how easily MMTK and BioPython could calculate the ϕ/ψ angles in radians. However, when drawing these plots it is conventional to use degrees:

```
import math
def degrees(rad_angle) :
    """Converts any angle in radians to degrees.

    If the input is None, then it returns None.
    For numerical input, the output is mapped to [-180,180]
    """
    if rad_angle is None :
        return None
    angle = rad_angle * 180 / math.pi
    while angle > 180 :
        angle = angle - 360
    while angle < -180 :
        angle = angle + 360
    return angle
```

## Downloads

You can download my python scripts to load a PDB file, calculate and save the Ramachandran angles (and the sample output files) here:

- Using MMTK:
  - [ramachandran_mmtk.py ![[Python Code\]](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/py.gif)](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/ramachandran_mmtk.py) - python script
  - [1HMP_mmtk.tsv](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/1HMP_mmtk.tsv) - sample output (plain text, tab separated variables)
- Using BioPython:
  - [ramachandran_biopython.py ![[Python Code\]](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/py.gif)](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/ramachandran_biopython.py) - python script
  - [1HMP_biopython.tsv](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/1HMP_biopython.tsv) - sample output (plain text, tab separated variables)

See also: [Parsing Lovell *et al.*'s Top 500 PDB files](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/top500/).

Next: Using these numbers to [draw the Ramachandran Plot with python](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/drawing/) (or [draw it with R](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/ramachandran/)).