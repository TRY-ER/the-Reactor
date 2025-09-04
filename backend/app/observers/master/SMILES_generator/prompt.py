QUERY = """
Following is the reaction text:

` {reaction_text} ` 
"""

RETRY_QUERY = """
You have made {invalid_count} invalid attempts.

Consider the errors and look into the syntax to come up with the right notation for SMILES representation.

Following is the reaction text:

` {reaction_text} ` 
"""

INSTRUCTIONS = """
# Details of SMILES notation

## SMILES - A Simplified Chemical Language
SMILES (Simplified Molecular Input Line Entry System) is a line notation (a typographical method using printable characters) for entering and representing molecules and reactions. Some examples are:


SMILES	| Name	
CC	| ethane	
[OH3+]	| hydronium ion
O=C=O	| carbon dioxide	
[2H]O[2H]	| deuterium oxide
C#N	| hydrogen cyanide	
[235U]	| uranium-235
CCN(CC)CC	| triethylamine	
F/C=C/F	| E-difluoroethene
CC(=O)O	| acetic acid	
F/C=C\F	| Z-difluoroethene
C1CCCCC1 | cyclohexane	
N[C@@H](C)C(=O)O  | L-alanine
c1ccccc1 | benzene	
N[C@H](C)C(=O)O	| D-alanine


## 3.2 SMILES Specification Rules
SMILES notation consists of a series of characters containing no spaces. Hydrogen atoms may be omitted (hydrogen-suppressed graphs) or included (hydrogen-complete graphs). Aromatic structures may be specified directly or in Kekulé form.
There are five generic SMILES encoding rules, corresponding to specification of atoms, bonds, branches, ring closures, and disconnections. Rules for specifying various kinds of isomerism are discussed in the following section, ISOMERIC SMILES.

### 3.2.1 Atoms
Atoms are represented by their atomic symbols: this is the only required use of letters in SMILES. Each non-hydrogen atom is specified independently by its atomic symbol enclosed in square brackets, [ ]. The second letter of two-character symbols must be entered in lower case. Elements in the "organic subset" B, C, N, O, P, S, F, Cl, Br, and I may be written without brackets if the number of attached hydrogens conforms to the lowest normal valence consistent with explicit bonds. "Lowest normal valences" are B (3), C (4), N (3,5), O (2), P (3,5), S (2,4,6), and 1 for the halogens. Atoms in aromatic rings are specified by lower case letters, e.g., aliphatic carbon is represented by the capital letter C, aromatic carbon by lower case c. Since attached hydrogens are implied in the absence of brackets, the following atomic symbols are valid SMILES notations.

C -> methane	(CH4)
P -> phosphine	(PH3)
N -> ammonia	(NH3)
S -> hydrogen sulfide	(H2S)
O -> water	(H2O)
Cl -> hydrochloric acid	(HCl)

Atoms with valences other than "normal" and elements not in the "organic subset" must be described in brackets.


[S] -> elemental sulfur
[Au] -> elemental gold

Within brackets, any attached hydrogens and formal charges must always be specified. The number of attached hydrogens is shown by the symbol H followed by an optional digit. Similarly, a formal charge is shown by one of the symbols + or -, followed by an optional digit. If unspecified, the number of attached hydrogens and charge are assumed to be zero for an atom inside brackets. Constructions of the form [Fe+++] are synonymous with the form [Fe+3]. Examples are:


[H+] -> proton
[Fe+2] -> iron (II) cation
[OH-] -> hydroxyl anion
[Fe++] -> iron (II) cation
[OH3+] -> hydronium cation
[NH4+] -> ammonium cation

### 3.2.2 Bonds
Single, double, triple, and aromatic bonds are represented by the symbols -, =, #, and :, respectively. Adjacent atoms are assumed to be connected to each other by a single or aromatic bond (single and aromatic bonds may always be omitted). Examples are:

CC	-> ethane -> (CH3CH3)
C=O	-> formaldehyde	-> (CH2O)
C=C	-> ethene	-> (CH2=CH2)
O=C=O	-> carbon dioxide	-> (CO2)
COC	-> dimethyl ether	-> (CH3OCH3)
C#N	-> hydrogen cyanide	-> (HCN)
CCO	-> ethanol	-> (CH3CH2OH)
[H][H]	-> molecular hydrogen	-> (H2)

For linear structures, SMILES notation corresponds to conventional diagrammatic notation except that hydrogens and single bonds are generally omitted. For example, 6-hydroxy-1,4-hexadiene can be represented by many equally valid SMILES, including the following three:


Structure -> Valid SMILES
 CH2=CH-CH2-CH=CH-CH2-OH -> [ C=CCC=CCO,
	C=C-C-C=C-C-O,
 	OCC=CCC=C]

### 3.2.3 Branches

Branches are specified by enclosing them in parentheses, and can be nested or stacked. In all cases, the implicit connection to a parenthesized expression (a "branch") is to the left. Examples are:

		
CCN(CC)CC	 -> Triethylamine	
Isobutyric acid -> CC(C)C(=O)O	
3-propyl-4-isopropyl-1-heptene -> C=CC(CCC)C(C(C)C)CCC


### 3.2.4 Cyclic Structures
Cyclic structures are represented by breaking one bond in each ring. The bonds are numbered in any order, designating ring opening (or ring closure) bonds by a digit immediately following the atomic symbol at each ring closure. This leaves a connected non-cyclic graph which is written as a non-cyclic structure using the three rules described above. Cyclohexane is a typical example:
C1CCCCC1

There are usually many different, but equally valid descriptions of the same structure, e.g., the following SMILES notations for 1-methyl-3-bromo-cyclohexene-1:

(a) CC1=CC(Br)CCC1

(b) CC1=CC(CCC1)Br


Many other notations may be written for the same structure, deriving from different ring closures. SMILES does not have a preferred entry on input; although (a) above may be simplest, others are just as valid.

A single atom may have more than one ring closure. This is illustrated by the structure of cubane in which two atoms have more than one ring closure:



Generation of SMILES for cubane: C12C3C4C1C5C4C3C25.

If desired, digits denoting ring closures can be reused. As an example, the digit 1 used twice in the specification:

O1CCCCC1N1CCCCC1
The ability to re-use ring closure digits makes it possible to specify structures with 10 or more rings. Structures that require more than 10 ring closures to be open at once are exceedingly rare. If necessary or desired, higher-numbered ring closures may be specified by prefacing a two-digit number with percent sign (%). For example, C2%13%24 is a carbon atom with a ring closures 2, 13, and 24 .

### 3.2.5 Disconnected Structures
Disconnected compounds are written as individual structures separated by a "." (period). The order in which ions or ligands are listed is arbitrary. There is no implied pairing of one charge with another, nor is it necessary to have a net zero charge. If desired, the SMILES of one ion may be imbedded within another as shown in the example of sodium phenoxide.

(a) [Na+].[O−]c1ccccc1

(b) c1cc([O−].[Na+])ccc1

Matching pairs of digits following atom specifications imply that the atoms are bonded to each other. The bond may be explicit (bond symbol and/or direction preceding the ring closure digit) or implicit (a nondirectional single or aromatic bond). This is true whether or not the bond ends up as part of a ring.

Adjacent atoms separated by dot (.) implies that the atoms are not bonded to each other. This is true whether or not the atoms are in the same connected component.

For example, C1.C1 specifies the same molecule as CC(ethane)

## 3.3 Isomeric SMILES
This section describes the SMILES rules used to specify isotopism, configuration about double bonds, and chirality. The term isomeric SMILES collectively refers to SMILES written using these rules.
The SMILES isomer specification rules allow chirality to be completely specified for any structure, if it is known. Unlike most existing chemical nomenclatures such as CIP and IUPAC, these rules are also designed to allow rigorous partial specification of chirality. Aside from use in macros, substructure searching, and other pattern matching operations, this is important because much of the world's available chemical information is known for structures with incompletely resolved chiralities (not all possible chiral centers are separated, known, or reported).

All isomer specification rules in SMILES are therefore optional. The absence of a specification for any attribute implies that the value of that attribute is unspecified.

### 3.3.1 Isotopic Specification
Isotopic specifications are indicated by preceding the atomic symbol with a number equal to the desired integral atomic mass. An atomic mass can only be specified inside brackets. For instance:

Smiles -> Name
[12C] -> carbon-12
[13C] -> carbon-13
[C]	-> carbon (unspecified mass)
[13CH4] -> C-13 methane

### 3.3.2 Configuration Around Double Bonds
Configuration around double bonds is specified by the characters / and \ which are "directional bonds" and can be thought of as kinds of single or aromatic (eg. default) bonds. These symbols indicate relative directionality between the connected atoms, and have meaning only when they occur on both atoms which are double bonded. For instance, the following SMILES are all valid for E- and Z-1,2-difluoroethene:

	
F/C=C/F	F/C=C\F

F\C=C\F	F\C=C/F

An important difference between SMILES chirality conventions and others such as CIP is that SMILES uses local chirality representation (as opposed to absolute chirality), which allows partial specifications. An example of this is illustrated below:
	
F/C=C/C=C/C	-> (completely specified)
F/C=C/C=CC -> (partially specified)
	
### 3.3.3. Configuration Around Tetrahedral Centers
SMILES uses a very general type of chirality specification based on local chirality. Instead of using a rule-based numbering scheme to order neighbor atoms of a chiral center, orientations are based on the order in which neighbors occur in the SMILES string. As with all other aspects of SMILES, any valid order is acceptable; the Daylight software is responsible for retaining the meaning of the chiral specification when the structure is modified or rearranged (e.g. to make the unique SMILES).
The simplest and most common kind of chirality is tetrahedral; four neighbor atoms are evenly arranged about a central atom, known as the "chiral center". If all four neighbors are different from each other in any way, mirror images of the structure will not be identical. The two mirror images are known as "enantiomers" and are the only two forms that a tetrahedral center can have. If two (or more) of the four neighbors are identical to each other, the central atom will not be chiral (its mirror images can be superimposed in space).

In SMILES, tetrahedral centers may be indicated by a simplified chiral specification (@ or @@) written as an atomic property following the atomic symbol of the chiral atom. If a chiral specification is not present for a chiral atom, its chirality is implicitly not specified. For instance:


(unspecified chirality)	-> [ NC(C)(F)C(=O)O , NC(F)(C)C(=O)O ]
(specified chirality) [ N[C@](C)(F)C(=O)O , N[C@@](F)(C)C(=O)O ]

Looking from the amino N to the chiral C (as the SMILES is written), the three other neighbors appear anticlockwise in the order that they are written in the top SMILES, N[C@](C)(F)C(=O)O (methyl-C, F, carboxy-C), and clockwise in the bottom one, N[C@@](F)(C)C(=O)O. The symbol "@" indicates that the following neighbors are listed anticlockwise (it is a "visual mnemonic" in that the symbol looks like an anticlockwise spiral around a central circle). "@@" indicates that the neighbors are listed clockwise (you guessed it, anti-anti-clockwise).

If the central carbon is not the very first atom in the SMILES and has an implicit hydrogen attached (it can have at most one and still be chiral), the implicit hydrogen is taken to be the first neighbor atom of the three neighbors that follow a tetrahedral specification. If the central carbon is first in the SMILES, the implicit hydrogen is taken to be the "from" atom. Hydrogens may always be written explicitly (as [H]) in which case they are treated like any other atom. In each case, the implied order is exactly as written in SMILES. Some of the valid SMILES for the alanine are:


Same representations below	
N[C@@]([H])(C)C(=O)O	
N[C@@H](C)C(=O)O	
N[C@H](C(=O)O)C	
[H][C@](N)(C)C(=O)O	
[C@H](N)(C)C(=O)O	

Same representations below
N[C@]([H])(C)C(=O)O
N[C@H](C)C(=O)O
N[C@@H](C(=O)O)C
[H][C@@](N)(C)C(=O)O
[C@@H](N)(C)C(=O)O

The chiral order of the ring closure bond is implied by the lexical order that the ring closure digit appears on the chiral atom (not in the lexical order of the "substituent" atom).


C[C@H]1CCCCO1
or
O1CCCC[C@@H]1C

### 3.3.4 General Chiral Specification
There are many kinds of chirality other than tetrahedral. The use of the "@" symbol described above is actually a special case of a general chiral specification syntax.
The general chiral specification used in SMILES has three parts: the @ symbol, followed by a two-letter chiral class indicator, followed by a numerical chiral permutation designator. A default chiral class is assigned to each degree (number of connections); the default class for four connections is tetrahedral (TH). Most chiralities have more than two possible choices; the choices are assigned from a table numerically. In most cases, the @1 designation means "anticlockwise around the axis represented by SMILES order" and @2 means "clockwise". Notations in the form "@@" and "@@@" are interpreted as "@2" and "@3" (analogous to "+++" meaning "+3"). The "@" and "@@" notations used above are shortcuts for the full specifications "@TH1" and "@TH2". In practice, full chiral specifications are not often needed.

SMILES handles the full range of chiral specification, including resolution of "reduced chirality" (where the number of enantiomers is reduced by symmetry) and "degenerate chirality" (where the center becomes non-chiral due to symmetrical substitution). As with other aspects of SMILES, the language guarantees the ability to specify exactly what is known, including partial specifications. The SMILES system will generate unique isomeric SMILES for any given specification, and substructure recognition will operate correctly on all types of chirality.

The rest of this section will be limited to discussing the following chiralities: tetrahedral, allene-like, square-planar, trigonal-bipyramidal, and octahedral. Although many more chiral classes can be handled by this system (it's table-driven), these five classes are very common in chemistry and cover most of the issues to be encountered in the remainder.

Tetrahedral. The tetrahedral class symbol is TH. This is the default chiral class for degree four. Possible values are 1 and 2. @TH1 (or just @) indicates that, looking from the first connected atom, the following three connected atoms are listed anticlockwise; @TH2 (or @@) indicates clockwise.

Allene-like. The allene-like class symbol is AL. This is the default chiral class for degree 2 (the chiral center is the central atom with two double bonds). Although substituted C=C=C structures are most common, C=C=C=C=C structures are also allene-like, as are any odd number of serially double-bonded atoms. Possible values are @AL1 (or just @) and @AL2 (or @@); these are interpreted by superimposing the substituted atoms and evaluating as per tetrahedral. Hydrogens attached to substituted allene-like atoms are taken to be immediately following that atom, as shown below:


	
OC(Cl)=[C@]=C(C)F -> OC(Cl)=[C@AL1]=C(C)F 
OC=[C@]=CF -> OC([H])=[C@AL1]=C([H])F

Square-planar. The square-planar class symbol is SP Possible values are @SP1, @SP2, and @SP3; this is not the default chiral class for degree four, so shorthand specifications are not allowed. Square-planar is also somewhat unusual in that the ideas of clockwise and anticlockwise do not apply.
	 
F[Po@SP1](Cl)(Br)I	(SP1 lists in a "U shape")
F[Po@SP2](Br)(Cl)I	(SP2 lists in a "4-shape")
F[Po@SP3](Cl)(I)Br	(SP3 lists in a "Z shape")

Trigonal-bipyramidal. The trigonal-bipyramidal class symbol is TB. This is the default chiral class for degree five. Possible values are @TB1 to @TB20. @TB1 (or just @) indicates that, when the SMILES is listed from one axial connection to the other, the three intermediate, equatorially-connected atoms are listed anticlockwise; @TB2 (or @@) indicates clockwise. This is illustrated below.

s[As@@](F)(Cl)(Br)C=O
O=C[As@](F)(Cl)(Br)S

Octahedral. The octahedral class symbol is OH. This is the default chiral class for degree six. Possible values are @OH1 to @OH30. @OH1 (or just @) indicates that, when the SMILES is listed from one axial connection to the other, the four intermediate, equatorially-connected atoms are listed anticlockwise; @OH2 (or @@) indicates clockwise. This is illustrated below.

S[Co@@](F)(Cl)(Br)(I)C=O
O=C[Co@](F)(Cl)(Br)(I)S

## 3.4 SMILES Conventions
Aside from the above rules, a small number of conventions are universally used in SMILES. These are briefly discussed below; for more detail, see the JCICS paper (ibid).

### 3.4.1 Hydrogens
Hydrogen atoms do not normally need to be specified when writing SMILES for most organic structures. The presence of hydrogens may be specified in three ways:

Implicitly.....for atoms specified without brackets, from normal valence assumptions.
Explicitly by count.....inside brackets, by the hydrogen count supplied; zero if unspecified.
As explicit atoms.....as [H] atoms.
There is no distinction between "organic" and "inorganic" SMILES nomenclature. One may specify the number of attached hydrogens for any atom in any SMILES. For example, propane may be entered as [CH3][CH2][CH3] instead of CCC.

There are four situations where specification of explicit hydrogen specification is required:

charged hydrogen, i.e. a proton, [H+];
hydrogens connected to other hydrogens, e.g., molecular hydrogen, [H][H];
hydrogens connected to other than one other atom, e.g., bridging hydrogens; and
isotopic hydrogen specifications, e.g. in heavy water, [2H]O[2H].

### 3.4.2 Aromaticity
Aromaticity must be deduced in a system such as SMILES which generates an unambiguous chemical nomenclature because of the fundamental requirement to characterize the symmetry of a molecule. Given effective aromaticity-detection algorithms, it is not necessary to enter any structure as aromatic if the user prefers to enter an aliphatic (Kekulé-like) structure. Entering structures as aromatic directly (i.e., by using lower case atomic symbols) provides a shortcut to accurate chemical specification and is closer to the mental molecular model used by most chemists.
The SMILES algorithm uses an extended version of Hueckel's rule to identify aromatic molecules and ions. To qualify as aromatic, all atoms in the ring must be sp2 hybridized and the number of available "excess" p-electrons must satisfy Hueckel's 4N+2 criterion. As an example, benzene is written c1ccccc1, but an entry of C1=CC=CC=C1 - cyclohexatriene, the Kekulé form - leads to detection of aromaticity and results in an internal structural conversion to aromatic representation. Conversely, entries of c1ccc1 and c1ccccccc1 will produce the correct anti-aromatic structures for cyclobutadiene and cyclooctatetraene, C1=CC=C1 and C1=CC=CC=CC=C1. In such cases the SMILES system looks for a structure that preserves the implied sp2 hybridization, the implied hydrogen count, and the specified formal charge, if any. Some inputs, however, may not only be incorrect but also impossible, such as c1cccc1. Here c1cccc1 cannot be converted to C1=CCC=C1 since one of the carbon atoms would be sp3 with two attached hydrogens. In such a structure alternating single and double bond assignments cannot be made. The SMILES system will flag this as an "impossible" input. Please note that only atoms on the following list can be considered aromatic: C, N, O, P, S, As, Se, and * (wildcard). In addition, exocyclic double bonds do not break aromaticity.

		
C1=COC=C1 -> c1cocc1
C1=CN=C[NH]C(=O)1 -> c1cnc[nH]c(=O)1	
C1=C*=CC=C1 -> c1c*ccc1
		
It is important to remember that the purpose of the SMILES aromaticity detection algorithm is for the purposes of chemical information representation only! To this end, rigorous rules are provided for determining the "aromaticity" of charged, heterocyclic, and electron-deficient ring systems. The "aromaticity" designation as used here is not intended to imply anything about the reactivity, magnetic resonance spectra, heat of formation, or odor of substances.

### 3.4.3 Aromatic Nitrogen Compounds
A short note is in order about aromatic nitrogens, a common source of confusion in chemical information systems. All three common types of aromatic nitrogen may be specified with the aromatic nitrogen symbol n. Archetypical examples are pyridine, pyridine-N-oxide, and pyrrole.

		
Pyridine -> n1ccccc1	

Pyridine-N-oxidei -> O=n1ccccc1   [O-][n+]1ccccc1	
Methyl and 1H-pyrrole -> Cn1cccc1    [nH]1cccc1

Note that the pyrrolyl nitrogen in 1H-pyrrole is written [nH] to distinguish this kind of nitrogen from a pyridyl-N. Alternative valid SMILES for 1H-pyrrole include [H]n1cccc1 (with explicit hydrogen) and N1C=CC=C1 (aliphatic form) all three input forms are equivalent.

### 3.4.4 Bonding Conventions

SMILES does not dictate which valence conventions should be used to model molecular structure. In fact, an advantage of using SMILES is its ability to describe various valence models of the same structure. Atoms may be connected and show charge separation as desired. For instance, nitromethane can be represented in SMILES as CN(=O)=O or as the charge separated C[N+](=O)[O-] (we tend to use the former for database work because it preserves symmetry). Both are "right" in the sense that they represent different, useful models of the substance. In general, when symmetry is not an issue, most chemists prefer charge-separated structures if they can avoid representing atoms in unusual valence states, e.g., diazomethane is written as C=[N+]=[N-] rather than C=[N]=[N].
Given one valence model of a structure, chemical database systems such as THOR and Merlin have the ability to retrieve data about that structure even if the data were stored under a different valence model of the structure. With such systems, the choice of valence conventions is not critical to either database design nor database query.

### 3.4.5 Tautomers
Tautomeric structures are explicitly specified in SMILES. There are no "tautomeric bond", "mobile hydrogen", nor "mobile charge" specifications. Selection of one or all tautomeric structures is left to the user and strongly depends on the application. Given one tautomeric form, most chemical information systems will report data for all known tautomers if needed. The role of SMILES is to specify exactly which tautomeric form is requested, and for which there are data. A simple example, with two possible tautomeric forms, is shown below:


2-pyridone -> O=c1[nH]cccc1         
2-pyridinol -> Oc1ncccc1
                
# Instruction

You are provided with set of texts that contain specifc reactions. Your job is to parse the reactant and products from it.
Once parsed convert that to corresponding valid SMILES and return us a ordered list of reactants and products SMILES.

"""