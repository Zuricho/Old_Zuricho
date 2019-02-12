# Calculating Ramachandran Angles

The **Biopython** library can be used to load PDB files and calculate the protein backbone's ϕ/ψ angles:

- Thomas Hamelryck's [Bio.PDB module in BioPython](https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ)
  Tolerant of odd PDB files, making it a safer bet.

Bio.PDB is tolerant of real world PDB files, it is the best choice for "data mining" tasks like calculating Ramachandran angles.



## ϕ/ψ using BioPython

[BioPython](http://www.biopython.org/) has a PDB parser, which includes its own vector classes which are used internally by the `Polypeptide` class to calculate a protein backbone's torsion angles. My first attempt looked like this:

```python
import Bio.PDB
for model in Bio.PDB.PDBParser().get_structure("1HMP", "1HMP.pdb") :
    for chain in model :
        poly = Bio.PDB.Polypeptide.Polypeptide(chain)
        print "Model %s Chain %s" % (str(model.id), str(chain.id)),
        print poly.get_phi_psi_list()
```

The less than perfect output, abbreviated:

```python
Model 0 Chain A
[(None, 2.3639878340439484), ..., (-1.4201969110624371, None)]
Model 0 Chain B
[(None, 2.1932190074065052), ..., (1.3986135567166249, None)]
Model 0 Chain  
[(None, None), ..., (None, None)]
```

This "worked" to an extent, but gave a extra polypeptide chain made out of the waters etc. The problem is the that code (above) ignored all the extra validation gained doing it using the `PPBuilder` or `CaPPBuilder` classes:

```python
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

```python
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

```python
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

```python
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



## Radians to Degrees

In the sections above we showed how easily BioPython could calculate the ϕ/ψ angles in radians. However, when drawing these plots it is conventional to use degrees:

```python
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



## Reference

1. https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/