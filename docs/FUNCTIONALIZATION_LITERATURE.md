<a id="funzionalizzazione-riferimenti-e-scelte-del-costruttore"></a>
# Functionalization: references and builder choices

Moro et al. (2024) describe GO modified with L-lysine, L-arginine and
L-methionine through epoxide opening under synthesis conditions at pH 9.
The XPS composition guides the model; functionalization also involves
partial reduction of the material. The study uses GAFF, TIP3P water and MD in
Amber 16; it reports RESP charges for glyphosate. These data support the
class of transformation, without automatically determining the reactive site,
protonation or charges of each new fragment in our generator.
Source: [article and official UNIBO record](https://cris.unibo.it/handle/11585/983516),
DOI [10.1016/j.mtchem.2024.101936](https://doi.org/10.1016/j.mtchem.2024.101936).
The main text was consulted; the DOCX supplement was not
retrieved during this verification.

Tunioli et al. (2023) describe amino-functionalized β-cyclodextrins with
C0/C2/C6/C12 linkers and their grafting onto GO through epoxide opening.
The control with unmodified βCD distinguishes the reactivity of the aminated
precursors under the conditions studied. The models are guided by XPS; the study
reports GAFF, charges obtained at the AM1 level, TIP3P, Amber 16 and
MM/GBSA analysis. Linker length and material composition are variables
in the study: they do not provide transferable parameters for every cyclodextrin or
arbitrary fragment.
Sources: [article](https://chemistry-europe.onlinelibrary.wiley.com/doi/full/10.1002/chem.202301854),
[official record and supplements](https://cris.unibo.it/handle/11585/948766).

<a id="regole-implementate"></a>
## Implemented rules

Our builder requires an explicitly selected neutral N/O/S site with a
transferable H. It breaks one epoxide C–O bond, forms the C–nucleophile bond and
transfers that H to the retained oxygen. Anti attack and a single attachment point
per opening are construction rules required for this project; they are not
presented as a reaction trajectory observed in the two articles.

For lysine, PEG, peptides or βCD, the atom is selected through the fragment map.
The α-N and ε-N of lysine remain distinct choices. Final protonation and the new
junction require their own parameters before MD; the program does not
automatically transfer intact-epoxide charges to the product.
