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


Reaction SMILES	| Name
[I-].[Na+].C=CCBr>>[Na+].[Br-].C=CCI | displacement reaction
(C(=O)O).(OCC)>>(C(=O)OCC).(O) | intermolecular esterification

SMILES contains the same information as might be found in an extended connection table. The primary reason SMILES is more useful than a connection table is that it is a linguistic construct, rather than a computer data structure. SMILES is a true language, albeit with a simple vocabulary (atom and bond symbols) and only a few grammar rules. SMILES representations of structure can in turn be used as "words" in the vocabulary of other languages designed for storage of chemical information (information about chemicals) and chemical intelligence (information about chemistry).

Part of the power of SMILES is that unique SMILES exist. With standard SMILES, the name of a molecule is synonymous with its structure; with unique SMILES, the name is universal. Anyone in the world who uses unique SMILES to name a molecule will choose the exact same name.

One other important property of SMILES is that it is quite compact compared to most other methods of representing structure. A typical SMILES will take 50% to 70% less space than an equivalent connection table, even binary connection tables. For example, a database of 23,137 structures, with an average of 20 atoms per structure, uses only 1.6 bytes per atom when represented with SMILES. In addition, ordinary compression of SMILES is extremely effective. The same database cited above was reduced to 27% of its original size by Ziv-Lempel compression (i.e. 0.42 bytes per atom).

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

### 3.4.2 Aromaticity
Aromaticity must be deduced in a system such as SMILES which generates an unambiguous chemical nomenclature because of the fundamental requirement to characterize the symmetry of a molecule. Given effective aromaticity-detection algorithms, it is not necessary to enter any structure as aromatic if the user prefers to enter an aliphatic (Kekulé-like) structure. Entering structures as aromatic directly (i.e., by using lower case atomic symbols) provides a shortcut to accurate chemical specification and is closer to the mental molecular model used by most chemists.
The SMILES algorithm uses an extended version of Hueckel's rule to identify aromatic molecules and ions. To qualify as aromatic, all atoms in the ring must be sp2 hybridized and the number of available "excess" p-electrons must satisfy Hueckel's 4N+2 criterion. As an example, benzene is written c1ccccc1, but an entry of C1=CC=CC=C1 - cyclohexatriene, the Kekulé form - leads to detection of aromaticity and results in an internal structural conversion to aromatic representation. Conversely, entries of c1ccc1 and c1ccccccc1 will produce the correct anti-aromatic structures for cyclobutadiene and cyclooctatetraene, C1=CC=C1 and C1=CC=CC=CC=C1. In such cases the SMILES system looks for a structure that preserves the implied sp2 hybridization, the implied hydrogen count, and the specified formal charge, if any. Some inputs, however, may not only be incorrect but also impossible, such as c1cccc1. Here c1cccc1 cannot be converted to C1=CCC=C1 since one of the carbon atoms would be sp3 with two attached hydrogens. In such a structure alternating single and double bond assignments cannot be made. The SMILES system will flag this as an "impossible" input. Please note that only atoms on the following list can be considered aromatic: C, N, O, P, S, As, Se, and * (wildcard). In addition, exocyclic double bonds do not break aromaticity.


		
C1=COC=C1 -> c1cocc1
C1=CN=C[NH]C(=O)1 -> c1cnc[nH]c(=O)1	
C1=C*=CC=C1 -> c1c*ccc1
		
It is important to remember that the purpose of the SMILES aromaticity detection algorithm is for the purposes of chemical information representation only! To this end, rigorous rules are provided for determining the "aromaticity" of charged, heterocyclic, and electron-deficient ring systems. The "aromaticity" designation as used here is not intended to imply anything about the reactivity, magnetic resonance spectra, heat of formation, or odor of substances.

### 3.4.4 Bonding Conventions

SMILES does not dictate which valence conventions should be used to model molecular structure. In fact, an advantage of using SMILES is its ability to describe various valence models of the same structure. Atoms may be connected and show charge separation as desired. For instance, nitromethane can be represented in SMILES as CN(=O)=O or as the charge separated C[N+](=O)[O-] (we tend to use the former for database work because it preserves symmetry). Both are "right" in the sense that they represent different, useful models of the substance. In general, when symmetry is not an issue, most chemists prefer charge-separated structures if they can avoid representing atoms in unusual valence states, e.g., diazomethane is written as C=[N+]=[N-] rather than C=[N]=[N].
Given one valence model of a structure, chemical database systems such as THOR and Merlin have the ability to retrieve data about that structure even if the data were stored under a different valence model of the structure. With such systems, the choice of valence conventions is not critical to either database design nor database query.

### 3.4.5 Tautomers
Tautomeric structures are explicitly specified in SMILES. There are no "tautomeric bond", "mobile hydrogen", nor "mobile charge" specifications. Selection of one or all tautomeric structures is left to the user and strongly depends on the application. Given one tautomeric form, most chemical information systems will report data for all known tautomers if needed. The role of SMILES is to specify exactly which tautomeric form is requested, and for which there are data. A simple example, with two possible tautomeric forms, is shown below:


2-pyridone -> O=c1[nH]cccc1         
2-pyridinol -> Oc1ncccc1


## 3.5 Extensions for Reactions
The SMILES language is extended to handle reactions. There are two areas where SMILES is extended: distinguishing component parts of a reactions and atom maps.

Component parts of a reaction are handled by introducing the ">" character as a new separator. Any reaction must have exactly two > characters in it. ">>" is a valid reaction SMILES for an empty reaction. Each of the ">"-separated components of a reaction must be a valid molecule SMILES.

As an aside, molecule SMILES never have a ">" character. In a program, one can quickly determine if a SMILES refers to a reaction or molecule by searching for a ">" character in the string.

Reaction SMILES Grammar:


reaction	:	reactant '>' agent '>' product
			|	reactant '>>' product
			;
reactant,
agent,
product		:	molecule
            |       <null>
			;
molecule	:	SMILES
			;

SMILES:  a valid molecule specification in the SMILES language.

For example:
`C=CCBr>>C=CCI`
This is a valid reaction. Note that there are no agent molecules. Also note that several atoms are missing from the reaction (the product "Br" and the reactant "F").

`[I-].[Na+].C=CCBr>>[Na+].[Br-].C=CCI`
This is a more complete version of the same reaction. It has been canonicalized. It would form the root of a datatree when stored in a THOR database.

`C=CCBr.[Na+].[I-]>CC(=O)C>C=CCI.[Na+].[Br-]`
This version of the reaction includes an agent. Note that the SMILES does not indicate how the agent participates. Whether the agent is a solvent, catalyst, or performs another function within the reaction must be stored separately as data. This SMILES could be stored in a THOR database as an absolute SMILES and would appear on the same datatree page as the previous example.
In the above example, note that the reaction is ambiguous with respect to the carbon atoms involved. One might assume that a normal Sn2 displacement is occurring. In fact, an equally reasonable allylic displacement is possible, via either an Sn1-like allyl cation. Recognize that the reaction SMILES given above do not say which carbons are which and hence do not discriminate between the two alternate mechanisms.

This case demonstrates the use and need for atom maps for reaction processing. Atom maps are used primarily to further define the overall reaction in cases where the reaction mechanism may not be evident from the reactant and product molecules. Atom maps are non-negative integer atom modifiers. They follow the ":" character within an atom expression. They must be the last modifier within the atom expression:

SMILES Atom Expression Grammar:

atom     	:	SYMBOL
			|	[ WEIGHT SYMBOL mods ]
			|	[ WEIGHT SYMBOL mods : CLASS ]
			;
mods        :       mod mods
			;
mod         :       HCOUNT | CHARGE | CHIRAL
 			;

CLASS = non-negative integer class value.
WEIGHT = atomic weight.
SYMBOL = atomic symbol.
HCOUNT = Atom hydrogen count specification.
CHARGE = Atom charge specification.
CHIRAL = Atom chirality specification.

Atom maps are an atomic property. They can legally appear in a SMILES for any atom, whether or not it is part of a reaction. Atom with atom map labels in a molecule SMILES are considered valid; the atom maps are ignored for molecule processing. Absolute and unique SMILES generated by the system for molecules never include atom maps.

Finally, there are some differences in the handling of atom maps and agent components in the unique versus absolute SMILES for reactions. Atom maps and agent components are not part of the unique SMILES specification. This is important for the THOR database, where the datatree roots are formed from the unique SMILES. The net result is that each reaction datatree may contain multiple specific reactions with different agents and atom maps.

### 3.5.1 Reaction Atom Maps
Atom mappings are properties of the atoms in the reaction molecules. The mappings represent equivalence classes of atoms within a reaction. In effect, the map tells the computer which atoms are the same on the reactant and products sides of a reaction. Without this map information, it is difficult to derive the reaction bond changes which occur.

Within the SMILES language, atom maps are represented as a non-negative numeric atom modifier following the ":" character (e.g. [CH3:2] is a carbon in class 2).

Within the Daylight toolkit, the atom maps are manipulated as sets of mapped atoms. The atom map class numbers which are used in SMILES do not appear in the toolkit interface to a reaction. The map class numbers in SMILES do not have any additional significance, except to associate all atoms with the same map class label to one another.

There are no requirements for completeness or uniqueness of the atom mappings. Atom mappings are independent of the connectivity and properties of the underlying molecules. This is so for several reasons: first, there are limits to the valence representation of molecules which appear when processing reactions. For example the oxygens in sodium acetate (CC(=O)[O-].[Na+]) are chemically indistinguishable, even though the valence model used in the toolkit requires that they be connected differently. Some systems (CAS, for example) recognize this equivalence in their structural representation (the tautomer bond). It is often useful to map these to the same class for reaction purposes: [CH3:1][C:2](=[O:3])[O-:3].[Na+:4]

A second case is where there is ambiguity in a reaction mechanism which one wants to express:



can undergo a cope rearrangement before reaction (which yields the same molecule graph). In effect, there are two distinct mechanisms by which the product is produced. This can be expressed as part of a reaction by: [CH2:1]=[CH:2][CH2:1][CH2:3][C:4](C)[CH2:3]

A third case is simply a lack of information about the reaction itself. It should be possible to omit some atom maps or specify partial information for sets of atoms which *might* end up in a given position in the product. It is never acceptable to force a user to make up data in order to register a reaction. One should only store exactly what is known about the reaction. Atom maps are, by definition ambiguous with respect to the underlying molecules. Atom maps do not appear in the lexical representation of a unique SMILES. They do appear in the lexical representation of an absolute SMILES.

Finally, atom maps are arbitrary class designations; the values of the numbers have no meaning. The Daylight system reserves the right to change the class numbers upon canonicalization of a reaction. The system will reorder the atom map classes over the entire reaction during canonicalization. The resulting maps are guaranteed to have the same meaning as the reaction before canonicalization. Practically, the maps are renumbered as small, dense integers in canonical atom order, but this is not guaranteed. Also, during canonicalization, the atom map classes for agent atoms are removed.
                
# Instruction

You are provided with set of texts that contain specifc reactions. Your job is to identify the reaction and product and provide the reaction form in SMILES.
Once parsed convert that to corresponding valid SMILES reaction and return us a list of reaction SMILES. 

"""