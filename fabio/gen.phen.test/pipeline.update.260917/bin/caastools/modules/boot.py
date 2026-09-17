#                      _              _     
#                     | |            | |    
#   ___ __ _  __ _ ___| |_ ___   ___ | |___ 
#  / __/ _` |/ _` / __| __/ _ \ / _ \| / __|
# | (_| (_| | (_| \__ \ || (_) | (_) | \__ \
#  \___\__,_|\__,_|___/\__\___/ \___/|_|___/


'''
A Convergent Amino Acid Substitution identification 
and analysis toolbox

Author:         Fabio Barteri (fabio.barteri@upf.edu)

Contributors:   Alejandro Valenzuela (alejandro.valenzuela@upf.edu)
                Xavier Farré (xfarrer@igtp.cat),
                David de Juan (david.juan@upf.edu).


MODULE NAME: boot.py
DESCRIPTION: bootstrap function
DEPENDENCIES: alimport.py, caas_id.py, pindex.py
CALLED BY: ct

'''


from modules.init_bootstrap import *
from modules.disco import process_position
from modules.caas_id import iscaas
from modules.alimport import *
from modules.hyper import calcpval_random
from modules.bootstrap_events import (
    EVENT_HEADER,
    infer_fixed_discovery_groups,
    natural_key,
    serialize_position_events,
    summarize_position_events,
)

from os.path import exists
import functools


# FUNCTION fetch_caas():
# fetches caas per each thing


# FUNCTION infer_resampled_group_sizes()
# Infers and validates the nominal FG/BG sizes from resampled traits.

def infer_resampled_group_sizes(resampled_traits):

    fg_traits = set(resampled_traits.trait2fg.keys())
    bg_traits = set(resampled_traits.trait2bg.keys())
    trait_names = fg_traits.union(bg_traits)

    if len(trait_names) == 0:
        raise ValueError("the resampled traits file contains no valid cycles")

    if fg_traits != bg_traits:
        raise ValueError("each resampling cycle must contain both foreground and background species")

    if len(trait_names) != resampled_traits.cycles:
        raise ValueError("resampling cycle names must be unique and every line must define a valid cycle")

    group_sizes = set()

    for trait in trait_names:
        fg_species = resampled_traits.trait2fg[trait]
        bg_species = resampled_traits.trait2bg[trait]

        if len(fg_species) == 0 or len(bg_species) == 0:
            raise ValueError("foreground and background groups cannot be empty")

        if len(fg_species) != len(set(fg_species)) or len(bg_species) != len(set(bg_species)):
            raise ValueError("a species cannot occur more than once in the same resampling group")

        if len(set(fg_species).intersection(set(bg_species))) > 0:
            raise ValueError("foreground and background groups must not overlap within a resampling cycle")

        group_sizes.add((len(fg_species), len(bg_species)))

    if len(group_sizes) != 1:
        formatted_sizes = ", ".join([str(x[0]) + "/" + str(x[1]) for x in sorted(group_sizes)])
        raise ValueError("all resampling cycles must have the same FG/BG sizes; found " + formatted_sizes)

    return list(group_sizes)[0]


# FUNCTION filter_positions_by_hypergeometric_pvalue()
# Excludes alignment positions whose random CAAS p-value exceeds the threshold.

def filter_positions_by_hypergeometric_pvalue(sliced_object, fg_size, bg_size, threshold):

    positions_before_filtering = len(sliced_object.d)
    retained_positions = []

    for position in sliced_object.d:
        pvalue = calcpval_random(position, sliced_object.genename, fg_size, bg_size)

        if pvalue <= threshold:
            retained_positions.append(position)

    sliced_object.d = retained_positions

    return positions_before_filtering, len(retained_positions)

# FUNCTION filter_for_gaps()
# filters a trait for its gaps

def filter_for_gaps(max_bg, max_fg, max_all, gfg, gbg):

    out = True

    all_g = gfg + gbg
    
    if max_all != "NO" and all_g > int(max_all):
        out = False

    elif max_fg != "NO" and gfg > int(max_fg):
        out = False

    elif max_bg != "NO" and gbg > int(max_bg):
        out = False

    return out


# FUNCTION filter_for_missing()
# filters a trait for its gaps

def filter_for_missings(max_m_bg, max_m_fg, max_m_all, mfg, mbg):

    out = True

    all_m = mfg + mbg
    
    if max_m_all != "NO" and all_m > int(max_m_all):
        out = False

    elif max_m_fg != "NO" and mfg > int(max_m_fg):
        out = False

    elif max_m_bg != "NO" and mbg > int(max_m_bg):
        out = False

    return out



def evaluate_bootstrap_hypothesis(
    trait,
    processed_position,
    maxgaps_fg,
    maxgaps_bg,
    maxgaps_all,
    maxmiss_fg,
    maxmiss_bg,
    maxmiss_all,
    minobserved_fg,
    minobserved_bg,
    admitted_patterns,
):
    """Return structured details for one positive hypothesis, else ``None``."""

    if trait not in processed_position.trait2aas_fg:
        return None
    if trait not in processed_position.trait2aas_bg:
        return None

    if len(processed_position.trait2ungapped_fg[trait]) < int(minobserved_fg):
        return None
    if len(processed_position.trait2ungapped_bg[trait]) < int(minobserved_bg):
        return None

    if maxgaps_fg != "NO" and processed_position.trait2gaps_fg[trait] > int(maxgaps_fg):
        return None
    if maxgaps_bg != "NO" and processed_position.trait2gaps_bg[trait] > int(maxgaps_bg):
        return None
    if maxgaps_all != "NO" and (
        processed_position.trait2gaps_fg[trait]
        + processed_position.trait2gaps_bg[trait]
        > int(maxgaps_all)
    ):
        return None

    if maxmiss_fg != "NO" and processed_position.trait2miss_fg[trait] > int(maxmiss_fg):
        return None
    if maxmiss_bg != "NO" and processed_position.trait2miss_bg[trait] > int(maxmiss_bg):
        return None
    if maxmiss_all != "NO" and (
        processed_position.trait2miss_fg[trait]
        + processed_position.trait2miss_bg[trait]
        > int(maxmiss_all)
    ):
        return None

    fg_amino_acids = sorted(set(processed_position.trait2aas_fg[trait]), key=natural_key)
    bg_amino_acids = sorted(set(processed_position.trait2aas_bg[trait]), key=natural_key)
    check = iscaas("/".join(["".join(fg_amino_acids), "".join(bg_amino_acids)]))
    if not check.caas or check.pattern not in admitted_patterns:
        return None

    species_to_amino_acid = {
        species: tagged_amino_acid.split("@", 1)[0]
        for species, tagged_amino_acid in processed_position.d.items()
    }
    fg_species_amino_acids = {
        species: species_to_amino_acid[species]
        for species in processed_position.trait2ungapped_fg[trait]
    }
    bg_species_amino_acids = {
        species: species_to_amino_acid[species]
        for species in processed_position.trait2ungapped_bg[trait]
    }

    return {
        "hypothesis": trait,
        "pattern": check.pattern,
        "fg_amino_acids": fg_amino_acids,
        "bg_amino_acids": bg_amino_acids,
        "fg_species_amino_acids": fg_species_amino_acids,
        "bg_species_amino_acids": bg_species_amino_acids,
    }


def caasboot(processed_position, genename, list_of_traits, maxgaps_fg, maxgaps_bg, maxgaps_all, maxmiss_fg, maxmiss_bg, maxmiss_all, minobserved_fg, minobserved_bg, cycles, admitted_patterns):

    valid_traits = sorted(
        set(list_of_traits)
        .intersection(processed_position.trait2aas_fg)
        .intersection(processed_position.trait2aas_bg),
        key=natural_key,
    )
    positive_hypotheses = []
    for trait in valid_traits:
        result = evaluate_bootstrap_hypothesis(
            trait,
            processed_position,
            maxgaps_fg,
            maxgaps_bg,
            maxgaps_all,
            maxmiss_fg,
            maxmiss_bg,
            maxmiss_all,
            minobserved_fg,
            minobserved_bg,
            admitted_patterns,
        )
        if result is not None:
            positive_hypotheses.append(result)

    output_traits = [result["hypothesis"] for result in positive_hypotheses]

    # Return the line
    
    position_name = genename + "@" + str(processed_position.position)
    count = str(len(output_traits))

    traitline = ",".join(output_traits)
    empval = str(int(count)/cycles)

    # The hypergeometric p-value describes the alignment position for the
    # nominal FG/BG group sizes. It is shared by all positive bootstrap cycles
    # reported in traitline and is printed immediately after that field.
    pvalue_string = str(processed_position.hypergeometric_pvalue)
    outline = "\t".join([position_name, count, str(cycles), empval, traitline, pvalue_string])

    return outline, positive_hypotheses
        

# FUNCTION disco_bootstrap()
# Launches the bootstrap in several lines. Returns a dictionary gene@position --> pvalue

def boot_on_single_alignment(trait_config_file, resampled_traits, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, min_fg_observed, min_bg_observed, the_admitted_patterns, output_file, event_output_file=None):


    # Step 3: processes the positions from imported alignment (process_position() from caas_id.py)
    processed_positions = list(map(functools.partial(process_position, multiconfig = resampled_traits, species_in_alignment = sliced_object.species), sliced_object.d))
    the_genename = sliced_object.genename
    print("caastools found", resampled_traits.cycles, "phenotype hypotheses")

    # Calculate the same positional hypergeometric p-value used by
    # --filter_significant, including when that prefilter is disabled. The
    # value is attached to the processed position so caasboot can serialize it.
    fg_size, bg_size = infer_resampled_group_sizes(resampled_traits)
    for raw_position, processed_position in zip(sliced_object.d, processed_positions):
        processed_position.hypergeometric_pvalue = calcpval_random(
            raw_position,
            the_genename,
            fg_size,
            bg_size,
        )

    # Step 4: extract the raw caas

    discovery_fg = None
    discovery_bg = None
    event_handle = None
    if event_output_file is not None:
        discovery_fg, discovery_bg = infer_fixed_discovery_groups(resampled_traits)
        event_handle = open(event_output_file, "w")
        print("\t".join(EVENT_HEADER), file=event_handle)

    with open(output_file, "w") as output_handle:
        for processed_position in processed_positions:
            line, positive_hypotheses = caasboot(
                processed_position=processed_position,
                list_of_traits=resampled_traits.alltraits,
                genename=the_genename,
                maxgaps_fg=max_fg_gaps,
                maxgaps_bg=max_bg_gaps,
                maxgaps_all=max_overall_gaps,
                maxmiss_fg=max_fg_miss,
                maxmiss_bg=max_bg_miss,
                maxmiss_all=max_overall_miss,
                minobserved_fg=min_fg_observed,
                minobserved_bg=min_bg_observed,
                admitted_patterns=the_admitted_patterns,
                cycles=resampled_traits.cycles,
            )
            print(line + "\t" + trait_config_file, file=output_handle)

            if event_handle is not None and positive_hypotheses:
                events = summarize_position_events(
                    processed_position,
                    positive_hypotheses,
                    discovery_fg,
                    discovery_bg,
                )
                rows = serialize_position_events(
                    gene=the_genename,
                    position=processed_position.position,
                    events=events,
                    positional_pvalue=processed_position.hypergeometric_pvalue,
                    trait_config=trait_config_file,
                )
                for row in rows:
                    print("\t".join(row), file=event_handle)

    if event_handle is not None:
        event_handle.close()

# FUNCTION pval()
# Returns a dictionary with the pvalue

def pval(bootstrap_result):
    with open(bootstrap_result) as h:
        thelist = h.read().splitlines()
    
    d = {}

    for line in thelist:
        try:
            c = line.split("\t")
            d[c[0]] = c[2]
        except:
            pass
    
    return d
