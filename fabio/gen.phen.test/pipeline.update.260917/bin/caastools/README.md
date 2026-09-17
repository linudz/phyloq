
This is derived from https://github.com/linudz/caastools

# CAAStools 1.0 - Software documentation.
# 1. Introduction to CAAStools

Amino acid substitutions that are consistent with phenotypic variation indicate that the gene product is potentially involved in the genetic determination of the trait. We define these cases as _Convergent Amino Acid Substitutions_ (CAAS).

It is possible to retrieve CAAS by scanning Multiple Sequence Alignments (MSA) of orthologous proteins or translated nucleotides. We will isolate those positions in which we can verify that species with diverging trait values converge to different amino acids. A very simple way to do this is to define two groups of species, isolate those positions in which they won’t share any amino acids, and test the statistical significance of this association.

Although the implementation of this strategy can be easily achieved through simple scripting, its scaling to proteome-size requires some additional effort in terms of code optimization. Also, CAAS analysis is usually validated through bootstrap-based approaches. All these operations together are computationally expensive, especially when brought to proteome-wide scale.

In recent years, our group worked on optimizing our in-house scripts for CAAS discovery and validation. CAAStools is the result of these efforts. A small suite of bioinformatics tools that allow the user to identify and validate CAAS analysis on MSA of orthologous proteins.


## 1.1 Suite overview 

CAAStools is a collection of 3 bioinformatics tools written in Python 3.9.4. **Figure 1** resumes the functioning of each tool. Globally, CAAStools relies on three main pieces of information that are required in different formats (see the input specification section of those paragraphs discussing each single tool). The **_discovery_** tool detects CAAS from a single amino acid MSA (protein or translated nucleotide). The **_resample_** tools elaborates virtual phenotype groups through different strategies (brownian motion simulation, random sorting of species or phylogeny-restricted resampling). The result of this tool can be submitted to the **_bootstrap_** tool that runs a bootstrap CAAS analysis on a single MSA.

[Link to **Figure 1**](https://drive.google.com/file/d/1sRTDl4ysuAriHcowLVC71zEQ31Bt3uHe/view)

The output of the discovery tool consists in a table reporting a list of CAAS associations. A corresponding p-value is calculated as the probability of randomly finding a CAAS association in that position (see our preprint for a full explanation about the statistical testing of CAAS).

## 1.2 Usage info and examples

The **caastools/examples** folder in the repository contains the files to run the examples described in this documentation. The supplementary information in our preprint ([Barteri et al., 2023](https://doi.org/10.1101/2022.12.14.520422)) describe the detection and the statistical testing of the gene BRCA2 from [Farré et al., 2021](https://doi.org/10.1093/molbev/msab219) analysis. The results and raw files about this example are in the folder **caastools/examples/longevity.2021**.

[Under construction]

# 2 Download, installation and updates.


## 2.1 Availability.

CAAStools is available as a GitHub repository

[https://github.com/linudz/caastools](https://github.com/linudz/caastools).  

The repository can be cloned through command line

`git clone github.com/linudz/caastools`

Also, you can clone it through the GitHub GUI client or download the code as a zipped file through your browser.


## 2.2 Installation.


The file _ct_ in the CAAStools Git folder is an executable python script that will launch the different tools of the suite. To run it system-wide, you can decide to either download the code to the usr/local/bin folder or to add the caastools folder to the $PATH variable.


### Updates.

Please, follow the github.com/linudz/caastools repository for code updates.


## 2.3 Dependencies and software requirements.

CAAStools is written in **Pyhton 3.9.4** and is compatible with any Python 3+ environment. Each tool has its own set of dependencies.


### Discovery and Bootstrap tools.

Biopython 1.79

Scipy (Biopython dependency)

Numpy (Biopython dependency)


### Resample tool

Dendropy 4+

Also, the _Brownian Motion_ _mode_ (`--mode bm`) for trait simulation relies on the simpermvec() function from R library RERconverge ([https://github.com/nclark-lab/RERconverge](https://github.com/nclark-lab/RERconverge)). This library will be requested only by the _resample_ tool when run in _Brownian Motion mode_ (`--mode bm`).

**Warning.** Apple M1 and M2 users might experience issues during RERconverge installation ([https://github.com/nclark-lab/RERconverge/issues/60](https://github.com/nclark-lab/RERconverge/issues/60)).


# 3. Discovery tool


## 3.1 Algorithm overview

CAAStools discovery identifies CAAS from an amino acid Multiple Sequence Alignment file (MSA). It is possible to view the inputs and the options through the `-h --help` command

`ct discovery -h/--help`

CAAStools discovery returns those positions in which amino acids differ between two groups of species. We call these groups **discovery groups**, distinguishing them into foreground (**FG**) and background (**BG**). The program will detect a CAAS in those MSA positions that will meet two conditions. First, the foreground and the background won’t share any amino acids. Second, all the species in at least one group need to share (or _converge_ to) the same amino acid. The two conditions define a set of 4 possible mutation **patterns**. _Table 3.1_ reports a set of possible mutation patterns, indicating which ones are accepted as CAAS, that are enumerated from 1 to 4.


<table>
  <tr>
   <td><strong>FG</strong>
   </td>
   <td><strong>BG</strong>
   </td>
   <td><strong>Difference between groups</strong>
   </td>
   <td><strong>Convergence within…</strong>
   </td>
   <td><strong>Pattern</strong>
   </td>
   <td><strong>Is it CAAS?</strong>
   </td>
  </tr>
  <tr>
   <td>K
   </td>
   <td>W
   </td>
   <td>YES
   </td>
   <td>Both groups
   </td>
   <td>Pattern 1
   </td>
   <td>YES
   </td>
  </tr>
  <tr>
   <td>K
   </td>
   <td>WE
   </td>
   <td>YES
   </td>
   <td>FG
   </td>
   <td>Pattern 2
   </td>
   <td>YES
   </td>
  </tr>
  <tr>
   <td>KE
   </td>
   <td>W
   </td>
   <td>YES
   </td>
   <td>BG
   </td>
   <td>Pattern 3
   </td>
   <td>YES
   </td>
  </tr>
  <tr>
   <td>KE
   </td>
   <td>WF
   </td>
   <td>YES
   </td>
   <td>None
   </td>
   <td>-
   </td>
   <td>NO
   </td>
  </tr>
  <tr>
   <td>K
   </td>
   <td>KE
   </td>
   <td>NO
   </td>
   <td>FG
   </td>
   <td>-
   </td>
   <td>NO
   </td>
  </tr>
  <tr>
   <td>KE
   </td>
   <td>K
   </td>
   <td>NO
   </td>
   <td>BG
   </td>
   <td>-
   </td>
   <td>NO
   </td>
  </tr>
</table>


**Table 3.1** - Mutation patterns.

Note that the **pattern 4** can be included as CAAS by user specification. Through the `--patterns` option, the user will be able to select the number of patterns to include in the output. By default, ct discovery returns the patterns 1,2 and 3.

The user can filter the result based on the maximum number of indels (or gaps, “-”) accepted per position or the maximum species missing in the alignment.

For each CAAS prediction, the program will calculate an empiric p-value that corresponds to the probability of finding the same set of mutational patterns with random species. This probability is calculated through the hypergeometric distribution.

For further details on the CAAStools discovery algorithm, please refer to our preprint on BioRXIV.


## 3.2 Formatting the inputs


### 3.2.1 The configuration file

`-t /--traitfile $config_file`

The **configuration file **or **config** is the file that we’ll use to tell the program which species we are comparing and how they are arranged into the FG and BG. It consists of a simple tab file containing the name of the species and a label indicating the corresponding group (**FG** = 1, **BG **= 0).

A config file is present in the examples/ folder (examples/conifig.tab)

Aotus_griseimembra	0

Avahi_peyrierasi	0

Callibella_humilis	0

Gorilla_beringei	1

Gorilla_gorilla	1

Macaca_thibetana	1

Mandrillus_leucophaeus	1

This config file will instruct the program to create two groups.


<table>
  <tr>
   <td>FG
   </td>
   <td>Gorilla_beringei, Gorilla_gorilla, Macaca_thibetana, Mandrillus_leucophaeus
   </td>
  </tr>
  <tr>
   <td>BG
   </td>
   <td>Aotus_griseimembra, Avahi_peyrierasi, Callibella_humilis
   </td>
  </tr>
</table>


Note that:



* The values are tab-separated and no further space is admitted
* The order of the species is irrelevant.


### 3.2.2 The amino acid MSA

The second fundamental input is the amino acid MSA file. CAAStools relies on Biopython 1.7 to import sequence files. Hence, the accepted formats are the ones specified in Biopython docs ( [https://biopython.org/wiki/AlignIO](https://biopython.org/wiki/AlignIO) ):

**clustal, emboss, fasta, fasta-m10, ig, maf, mauve, msf, nexus, phylip, phylip-sequential, phylip-relaxed, stockholm**

By default, ct discovery will read the MSA as a clustal file. To specify a different format, we’ll need to specify it through the `--fmt` option (e.g. `--fmt phylip-relaxed`). In the examples folder, the examples/MSA directory contains an MSA in different formats.


### 3.2.3 A note on name consistency

CAAStools will associate the sequence in the MSA to the species by name identity. The program will save the name of each species in FG and BG groups in string variables, and will select those sequences in the MSA whose ID field will coincide with one of the species in the config file. **Please, format your config file and MSA in order to match the sequence IDs with the name of the species in the config file**.


## 3.3 Gaps and missing species filtering 

CAAS results can be filtered for a maximum of gaps or missing species. In this, we define as a “gap” the presence of an indel which is indicated with the “-” character. A missing species will be a species that is mentioned in the config file but it is not found in the MSA. This situation can occur when we iterate the analysis over different MSAs with variable coverage.

By default, CAAStools discovery accepts a maximum of n-1 gaps or missing species per group, where n is the size of the group. This means that the program will need the presence of at least one species per group to verify the conditions for CAAS assignment. The user can decide to limit the number of gaps and missing species per group, or to skip those positions in which gaps represent more than a maximum percentage of total symbols (default=50%). The following tables report the different options for gaps and missing species.


### 3.3.1 Filtering for gaps


<table>
  <tr>
   <td>Max background gaps
   </td>
   <td>--max_bg_gaps
   </td>
   <td>Filter by number of gaps in the background.
   </td>
   <td>Default: No filter
   </td>
  </tr>
  <tr>
   <td>Max foreground gaps
   </td>
   <td>--max_fg_gaps
   </td>
   <td>Filter by number of gaps in the foreground.
   </td>
   <td>Default: No filter
   </td>
  </tr>
  <tr>
   <td>Max overall gaps
   </td>
   <td>--max_gaps
   </td>
   <td>Filter by total number of gaps
   </td>
   <td>Default: No filter
   </td>
  </tr>
  <tr>
   <td>Max gaps per position
   </td>
   <td>--max_gaps_per_position
   </td>
   <td>Filter by number of gaps per position.
   </td>
   <td>Default: 0.5 (50%)
   </td>
  </tr>
</table>



### 3.3.2 Filtering for missing species


<table>
  <tr>
   <td>Max background missing species
   </td>
   <td>--max_bg_miss
   </td>
   <td>Filter by number of missing species in the background.
   </td>
   <td>Default: No filter
   </td>
  </tr>
  <tr>
   <td>Max foreground missing species
   </td>
   <td>--max_fg_miss
   </td>
   <td>Filter by number of missing species in the foreground.
   </td>
   <td>Default: No filter
   </td>
  </tr>
  <tr>
   <td>Max overall missing species
   </td>
   <td>--max_miss
   </td>
   <td>Filter by total number of missing species
   </td>
   <td>Default: No filter
   </td>
  </tr>
</table>

## 3.4 The output

CAAStools discovery will output a tab file with all the CAAS found in one single MSA.


<table>
  <tr>
   <td><strong>Column</strong>
   </td>
   <td><strong>Header</strong>
   </td>
   <td><strong>Description</strong>
   </td>
  </tr>
  <tr>
   <td>1
   </td>
   <td>Gene
   </td>
   <td>The name of the gene (from MSA filename)
   </td>
  </tr>
  <tr>
   <td>2
   </td>
   <td>Trait
   </td>
   <td>The name of the trait (from binary config file)
   </td>
  </tr>
  <tr>
   <td>3
   </td>
   <td>Position
   </td>
   <td>The position in the MSA (0-based)
   </td>
  </tr>
  <tr>
   <td>4
   </td>
   <td>Substitution
   </td>
   <td>The substitution FG/BG
   </td>
  </tr>
  <tr>
   <td>5
   </td>
   <td>Pvalue
   </td>
   <td>The p-value from hypergeometric distribution
   </td>
  </tr>
  <tr>
   <td>6
   </td>
   <td>Scenario
   </td>
   <td>The mutational pattern (see “patterns” table in<em> 3.1 - Algorithm overview</em>)
   </td>
  </tr>
  <tr>
   <td>7
   </td>
   <td>FFGN
   </td>
   <td>Species <strong>F</strong>ound in <strong>FG</strong>: <strong>N</strong>umber. Number of species found in the FG group (it excludes those ones having an indel). 
   </td>
  </tr>
  <tr>
   <td>8
   </td>
   <td>FBGN
   </td>
   <td>Species <strong>F</strong>ound in <strong>BG</strong>: <strong>N</strong>umber. Number of species found in the BG group (it excludes those ones having an indel). 
   </td>
  </tr>
  <tr>
   <td>9
   </td>
   <td>GFG
   </td>
   <td>Number of <strong>G</strong>aps in the <strong>FG</strong>.
   </td>
  </tr>
  <tr>
   <td>10
   </td>
   <td>GBG
   </td>
   <td>Number of <strong>G</strong>aps in the <strong>BG</strong>.
   </td>
  </tr>
  <tr>
   <td>11
   </td>
   <td>MFG
   </td>
   <td>Number of <strong>M</strong>issing species in the <strong>FG</strong>.
   </td>
  </tr>
  <tr>
   <td>12
   </td>
   <td>MBG
   </td>
   <td>Number of <strong>M</strong>issing species in the <strong>BG</strong>.
   </td>
  </tr>
  <tr>
   <td>13
   </td>
   <td>FFG
   </td>
   <td>List of species <strong>F</strong>ound in the <strong>FG</strong>. Comma-separated.
   </td>
  </tr>
  <tr>
   <td>14
   </td>
   <td>FBG
   </td>
   <td>List of species <strong>F</strong>ound in the <strong>FG</strong>. Comma-separated.
   </td>
  </tr>
  <tr>
   <td>15
   </td>
   <td>MS
   </td>
   <td>List of missing species. Comma-separated.
   </td>
  </tr>
</table>



## 


## 3.5 Examples

Run ct discovery with example alignment (`phylip-relaxed` format, that has to be specified).

`ct discovery -a examples/MSA/primates.msa.pr -t examples/config.tab -o examples/discovery.output.usr.example --fmt phylip-relaxed`


# 4 Resample tool


## 4.1 Algorithm Overview

CAAStools resample (ct resample) elaborates a set of resampled discovery groups for bootstrap analysis. The global options for the tool can be fetched through the help command:

`ct resample -h/--help`

The simulation of discovery groups is propaedeutic to bootstrap analysis

can consist in a simple randomization, a randomization that is restricted to some parts of the phylogeny, or be based on brownian-motion trait evolution simulation. The user will indicate one of these three strategies, the size of the resampled FG and BG groups and the number of simulation cycles. The program outputs a tab file in 


### 4.1.1 Random simulation strategy

`ct resample --mode random`

In this case, the simulation will consist in the bare random sorting of species into a pair of FG/BG discovery groups. 


### 4.1.2 Phylogeny-restricted random simulation strategy.

`ct resample --mode random --limit_by_group $groupfile`

In this case, the simulation is based on the random choice of species, but is limited to the families that are present in a config file provided as a template. A further file, the species file, specifies the composition of the families. The random scooping takes into account the number of groups (or families) present in the template groups and will replicate that composition. For instance, if our template FG group consists of 3 species from group A and 2 species from groupB, the randomisation will follow this pattern. In each cycle, the program scoops 3 random species from group A and 2 random species from group B.


### 4.1.3 Brownian motion based simulation strategy.

`ct resample --mode bm`

This strategy resamples neutral evolution by brownian motion simulation. The FG/BG group size is defined by a template config file. The R function

`simpermvec()`

from R library RERconverge ([https://github.com/nclark-lab/RERconverge](https://github.com/nclark-lab/RERconverge)) will perform a brownian motion simulation of neutral evolution. FG and BG will be defined as the n-th species with higher values and the m-th species with lower values, where n and m are the size of FG and BG respectively.


## 4.2 Inputs per simulation strategy

Each simulation strategy will require a specific set of input files. The following table reports all the inputs that are needed for each strategy.


<table>
  <tr>
   <td><strong>Input File</strong>
   </td>
   <td><strong>Random</strong>
   </td>
   <td><strong>Random (Phylogeny restricted)</strong>
   </td>
   <td><strong>Brownian Motion</strong>
   </td>
  </tr>
  <tr>
   <td>Phylogenetic tree in newick format.
<p>
-p/--phylogeny
   </td>
   <td>Mandatory
   </td>
   <td>Mandatory
   </td>
   <td>Mandatory
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>Config File as a template
<p>
--bytemp
   </td>
   <td>Can be replaced by -f/--fg_size and -b/--bg_size for FG/BG size definition
   </td>
   <td>Mandatory
   </td>
   <td>Mandatory
   </td>
  </tr>
  <tr>
   <td>Trait values
<p>
--traitvalues
   </td>
   <td>NO
   </td>
   <td>NO
   </td>
   <td>Mandatory. It is used by the program to shuffle the
   </td>
  </tr>
</table>

The user will need to provide a phylogenetic tree in Newick format ([https://evolution.genetics.washington.edu/phylip/newicktree.html](https://evolution.genetics.washington.edu/phylip/newicktree.html)) to tell the program on which species base its simulation. The binary configuration file is required for the phylogeny restricted and Brownian motion strategies. The phylogeny restricted strategy limits trait randomisation to some specific families.

The resample tool outputs 1000 resampled traits by default. The user can decide the number of cycles through the --cycles option:

`--cycles 10`

`--cycles 100`

`--cycles 1000`


## 4.3 The output

The output consists in a tab file with three columns.

**Column 1**: The name of the cycle, indicated in the b_numberofcycle format

**Column 2**: The FG species (comma-separated)

**Column 3**: The BG species (comma-separated)

Here’s an example of the first ten lines of a resampled traits output file.

  b_1	Cercopithecus_mitis,Alouatta_palliata,Pan_troglodytes,Colobus_polykomos	Saimiri_sciureus,Alouatta_puruensis,Trachypithecus_crepusculus,Cercocebus_torquatus

  b_2	Macaca_fuscata,Papio_cynocephalus,Cercopithecus_petaurista,Cheracebus_lucifer	Trachypithecus_phayrei,Tarsius_wallacei,Lophocebus_aterrimus,Macaca_silenus

  b_3	Saimiri_oerstedii,Nomascus_concolor,Lemur_catta,Saguinus_oedipus	Pongo_abelii,Indri_indri,Eulemur_albifrons,Eulemur_fulvus

  b_4	Saimiri_ustus,Eulemur_rubriventer,Leontocebus_nigricollis,Macaca_mulatta	Semnopithecus_hypoleucos,Mus_musculus,Otolemur_garnettii,Eulemur_macaco

  b_5	Cacajao_hosomi,Alouatta_puruensis,Saimiri_macrodon,Pithecia_mittermeieri	Ateles_belzebuth,Macaca_maura,Prolemur_simus,Trachypithecus_laotum

  b_6	Eulemur_coronatus,Trachypithecus_geei,Hapalemur_griseus,Prolemur_simus	Eulemur_flavifrons,Macaca_fuscata,Trachypithecus_pileatus,Cercopithecus_mona

  b_7	Cheracebus_lugens,Cercocebus_lunulatus,Hapalemur_occidentalis,Ateles_chamek	Tarsius_lariang,Cercopithecus_neglectus,Cercopithecus_diana,Macaca_thibetana

  b_8	Perodicticus_potto,Pan_paniscus,Cacajao_hosomi,Lepilemur_ankaranensis	Cercocebus_chrysogaster,Papio_anubis,Callimico_goeldii,Plecturocebus_miltoni

  b_9	Macaca_leonina,Cercopithecus_diana,Propithecus_diadema,Macaca_siberu	Galagoides_demidovii,Alouatta_seniculus,Propithecus_coquereli,Plecturocebus_dubius

  b_10	Leontopithecus_rosalia,Papio_cynocephalus,Hylobates_agilis,Mus_musculus	Cheracebus_regulus,Cercopithecus_mitis,Lophocebus_aterrimus,Ateles_marginatus


## 4.4 Examples {#4-4-examples}


### Ex.1 Resampling based on random selection of species {#ex-1-resampling-based-on-random-selection-of-species}

**With input fg/bg size (-f and -b options)**

 `ct resample -p examples/phylogeny.nw -f 5 -b 4 -m random --cycles 500 -o test/resample/random.resampling.tab`

**By template (binary config)**

`ct resample -p examples/phylogeny.nw --bytemp examples/config.tab -m random --cycles 500 -o test/resample/random.resampling.bytemplate.tab`

**Phylogeny restricted (must go by template)**

`ct resample -p examples/phylogeny.nw --bytemp examples/config.tab -m random --limit_by_group test/sp2fam.210727.tab --cycles 500 -o test/resample/random.resampling.bytemplate.tab`


### Ex.2 resampling based on BM 

**Template and trait values mandatory**

`ct resample -p examples/phylogeny.nw --bytemp examples/config.tab -m random --cycles 500  --traitvalues examples/traitvalues.tab -o test/resample/BM.resampling.tab`




# 5 Bootstrap Tool

**NOTE**: The comparison of the different statistical testings are described in Supplementary Material 3 in [CAAStools manuscript](https://doi.org/10.1101/2022.12.14.520422).

## 5.1 Algorithm overview 

The bootstrap tool is designed to repeat the CAAS discovery on a large number of discovery groups. In this case, the discovery groups are defined through the output file of the resample tool, in which each line represents a single cycle (see previous paragraph). The program will scan an MSA and will count the number of cycles that return a CAAS in that position.


## 5.2 The inputs

`-s $resampled_trait (output of ct resample)`

`-a $MSA`

### 5.2.1 Optional hypergeometric significance filter

By default, the bootstrap tool evaluates every alignment position retained by
the standard alignment and gap filters. The `--filter_significant` option can
be used to exclude positions whose hypergeometric CAAS p-value is greater than
a user-defined threshold:

`--filter_significant 0.05`

The foreground and background sizes are inferred from the resampled traits
file. All resampling cycles must contain the same non-zero FG/BG sizes, unique
species within each group, and no overlap between foreground and background.
The p-value is calculated once per alignment position using those nominal
group sizes, before the bootstrap cycles are evaluated.

The default value is the string `no`, which preserves the unfiltered bootstrap
behaviour. The aliases `--filter-significant` and `-fs` are also accepted.

### 5.2.2 Minimum observed foreground and background species

Bootstrap cycles can be required to have a minimum number of usable species
in each group at every evaluated alignment position:

`--min_fg_observed 3 --min_bg_observed 3`

An observed species is a member of the resampled FG or BG group that is
present in the alignment and non-gapped at that position. If either group is
below its requested minimum, that cycle is not counted as a CAAS at the
position. The position remains in the bootstrap output and may still be
positive for other cycles with sufficient coverage.

Both options default to `1` for backwards compatibility. Their hyphenated
aliases, `--min-fg-observed` and `--min-bg-observed`, are also accepted. A
minimum must be a positive integer and cannot exceed the corresponding group
size inferred from the resampled traits file.

### 5.2.3 Species- and amino-acid-aware event summaries

The historical bootstrap count treats every positive resampling row as one
positive cycle. This is useful for discovery, but the resampling rows may
overlap strongly and should not automatically be interpreted as independent
biological observations. The optional event summary retains the positive
resampling identifiers as **hypotheses** and collapses compatible hypotheses
into species-level amino-acid events:

`--summarize_species yes`

The companion table is written to `OUTPUT.events.tsv`. An explicit path can be
provided with `--event_output FILE`; providing that option also enables the
summary. Hyphenated aliases (`--summarize-species` and `--event-output`) are
accepted.

An event signature contains the set of FG amino acids and the disjoint set of
BG amino acids observed among its positive hypotheses. Two signatures are
compatible when the union of their FG amino acids remains disjoint from the
union of their BG amino acids. Compatible signatures are combined using all
maximal compatible groups, rather than a greedy and order-dependent merge.
Incompatible signatures at the same position remain separate events. Events
are ranked by balanced FG/BG discovery support, total unique-species support,
nominal event p-value, and signature; the first is marked as primary.

Each species is counted once per event even if it occurs in many positive
hypotheses. Counts are reported against two denominators:

- **observed species**: discovery-pool species present in the alignment with a
  non-gap amino acid at that position;
- **discovery species**: the complete fixed FG or BG pool inferred from all
  hypotheses in the resampling file.

The event summary requires fixed sides: a species may occur in many
hypotheses, but it cannot be FG in one hypothesis and BG in another. This
constraint is checked only when event output is requested, so ordinary random
bootstrap analyses remain backwards compatible.

For each event, CAAStools also reports species whose amino acid conflicts with
the event signature (a BG species carrying an event FG amino acid or an FG
species carrying an event BG amino acid). The nominal event p-value is an
exact conditional label test. It preserves the observed FG/BG sample sizes
and evaluates the minimum of the two oriented signature-match fractions. Its
calculation is reduced to three sufficient amino-acid categories, avoiding an
expensive enumeration of individual permutations. Because the event signature
was selected from the same data, this p-value is explicitly nominal;
proteome-wide correction and alignment-cluster filtering remain downstream.


## 5.3 Output

The output consists of a headerless, tab-separated file with seven columns.

Column 1: Position

Column 2: Number of resamples returning a CAAS in the position

Column 3: Number of cycles

Column 4: Bootstrap value

Column 5: Cycles with positive CAAS

Column 6: Hypergeometric p-value for the alignment position and the nominal
foreground/background group sizes

Column 7: Trait-template filename

The p-value in column 6 belongs to the position, so one value follows the
complete comma-separated cycle field even when several cycles are positive.

When event summaries are enabled, a second headered, tab-separated long table
contains one row per position/event. Its fields include:

- event identifier and primary-event flag;
- merged FG and BG amino-acid signature;
- dominant FG and BG amino acid and its proportion among the unique species
  supporting that side of the event (for example `A=0.7`); all co-dominant
  residues are retained when their counts are tied;
- positive hypothesis identifiers and the CAAS pattern of the final merged
  event (1 for one-vs-one, 2 for one-vs-many, 3 for many-vs-one, and 4 for
  many-vs-many amino acids);
- unique supporting species and their amino acids for FG and BG;
- amino-acid counts based on unique species, never hypothesis occurrences;
- observed and complete discovery species, counts, and denominators;
- support fractions over observed and discovery species;
- balanced FG/BG support and total unique-species support;
- species and amino acids conflicting with the event signature;
- nominal event p-value, positional hypergeometric p-value, and trait template.

Only positive positions appear in the event table. The original bootstrap
table still contains all positions that reached the bootstrap stage, including
rows with zero positive hypotheses.


## 5.4 Examples


### Bootstrap from random resampled traits.

`ct bootstrap -s test/resample/random.resampling.tab -t examples/config.tab -a examples/MSA/primates.msa.pr -o examples/random.bootstrap.tab --fmt phylip-relaxed`

### Bootstrap restricted to hypergeometrically significant positions.

`ct bootstrap -s test/resample/random.resampling.tab -t examples/config.tab -a examples/MSA/primates.msa.pr -o examples/random.bootstrap.significant.tab --fmt phylip-relaxed --filter_significant 0.05`

### Bootstrap requiring at least three observed species per group.

`ct bootstrap -s test/resample/random.resampling.tab -t examples/config.tab -a examples/MSA/primates.msa.pr -o examples/random.bootstrap.minimum3.tab --fmt phylip-relaxed --min_fg_observed 3 --min_bg_observed 3`

### Bootstrap with a species/amino-acid event companion table.

`ct bootstrap -s fixed-sides.resampling.tsv -t examples/config.tab -a examples/MSA/primates.msa.pr -o examples/fixed.bootstrap.tsv --fmt phylip-relaxed --min_fg_observed 3 --min_bg_observed 3 --summarize_species yes`

## 5.5 Pooled discovery

`pooled-discovery` is the species-aware successor to the fixed-side bootstrap
workflow. Instead of requiring an externally generated resampling table, it
accepts one complete discovery-pool config through `-t`. The input contains
one species and one binary side per line:

```text
species_A	1
species_B	1
species_C	0
species_D	0
```

CAAStools validates unique species and non-empty FG/BG pools, enumerates all
unique combinations of the requested FG and BG subset sizes, selects the
requested number of paired comparisons, saves those realized hypotheses, and
runs the existing position scanner and event summarizer in one command.

The principal options are:

- `--fg_size` / `--fg-size`: number of FG species per hypothesis (default 4);
- `--bg_size` / `--bg-size`: number of BG species per hypothesis (default 4);
- `--comparisons max|N`: number of unique paired hypotheses. The default
  `max` uses every possible comparison;
- `--seed`: seed controlling a smaller deterministic selection (default
  260811);
- `--hypotheses_output auto|none|FILE`: where to save the realized headerless
  resampling-format table. `auto`, the default, creates
  `OUTPUT.hypotheses.tsv`;
- `--event_output FILE`: optional explicit event-table location. Species event
  summarization is enabled by default for `pooled-discovery`.

For seven FG and six BG species with four selected from each side, there are:

```text
choose(7, 4) * choose(6, 4) = 35 * 15 = 525
```

Thus `--comparisons max` generates 525 hypotheses, while
`--comparisons 100` generates 100 unique hypotheses without replacement. For
a selection smaller than the maximum, candidates are ranked by a SHA-256 key
constructed from the seed and the two species subsets. This is deterministic
across runs and independent of Python's `random.sample` implementation. It
does not reproduce an older R `sample()` selection merely because the numeric
seed is the same; the realized hypotheses file is therefore always saved for
provenance.

The complete pools, not only the species occurring in a selected subset, are
retained as the FG/BG discovery denominators in the event table. The alignment
slicer still uses the per-hypothesis subset sizes, making its candidate-position
behavior equivalent to an external four-vs-four resampling table.

Example using 100 of the 525 possible longevity comparisons:

```bash
ct pooled-discovery \
  -a GENE.phy \
  -t longevity.full-pools.caas.cfg \
  -o GENE.pooled100.caas.tsv \
  --event-output GENE.pooled100.caas.events.tsv \
  --hypotheses-output longevity.pooled100.hypotheses.tsv \
  --fmt phylip-relaxed \
  --fg-size 4 --bg-size 4 \
  --comparisons 100 --seed 260811 \
  --min-fg-observed 3 --min-bg-observed 3
```

The original `bootstrap` command remains available and unchanged for analyses
that already have a resampling table or allow species to change sides.

5. License

This software is licensed under GNU General Public License. The kind of license is to be decided with UPF.


# 6. How to cite

Barteri, F., Valenzuela, A., Farré, X., de Juan, D., Muntané, G., Esteve-Altava, B., & Navarro, A. (2023). CAAStools: a toolbox to identify and test Convergent Amino Acid Substitutions. Bioinformatics, 39(10), btad623. 
[Open Access](https://doi.org/10.1093/bioinformatics/btad623)

# 7. Questions and troubleshooting 

You can ask your questions through the [discussions section](https://github.com/linudz/caastools/discussions) of CAAStools github. Also, you can contact Fabio Barteri at Pompeu Fabra University / BBRC [fabio.barteri@upf.edu](mailto:fabio.barteri@upf.edu)
