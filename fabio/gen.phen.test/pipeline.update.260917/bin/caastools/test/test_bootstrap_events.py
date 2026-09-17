#!/usr/bin/env python3

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


CAASTOOLS_DIR = Path(__file__).resolve().parents[1]
CT = CAASTOOLS_DIR / "ct"

sys.path.insert(0, str(CAASTOOLS_DIR))

from modules.bootstrap_events import (
    EVENT_HEADER,
    classify_event_pattern,
    exact_signature_separation_pvalue,
    format_dominant_amino_acid,
    summarize_position_events,
)


class BootstrapEventTests(unittest.TestCase):

    def test_final_event_pattern_uses_merged_amino_acid_signature(self):
        self.assertEqual(classify_event_pattern({"A"}, {"V"}), "1")
        self.assertEqual(classify_event_pattern({"A"}, {"I", "V"}), "2")
        self.assertEqual(classify_event_pattern({"A", "G"}, {"V"}), "3")
        self.assertEqual(classify_event_pattern({"A", "G"}, {"I", "V"}), "4")

    def test_dominant_amino_acid_uses_unique_species_and_preserves_ties(self):
        self.assertEqual(
            format_dominant_amino_acid({"species_1": "A", "species_2": "A", "species_3": "G"}),
            "A=0.666666666667",
        )
        self.assertEqual(
            format_dominant_amino_acid({"species_1": "A", "species_2": "G"}),
            "A=0.5,G=0.5",
        )

    def write_compatible_inputs(self, directory):
        directory = Path(directory)
        alignment = directory / "compatible.phy"
        trait_config = directory / "trait.cfg"
        resampled_traits = directory / "resampled.tsv"

        # The first three hypotheses are positive and have signatures A/V,
        # A/{I,V}, and {A,G}/V. They can be merged into one {A,G}/{I,V}
        # event. The fourth hypothesis completes the seven-FG/six-BG fixed
        # discovery pools but is pattern 4 and therefore not positive under
        # the default admitted patterns.
        alignment.write_text(
            "13 1\n"
            "Aotus_trivirgatus         A\n"
            "Ateles_geoffroyi          A\n"
            "Cercocebus_torquatus      A\n"
            "Eulemur_fulvus            A\n"
            "Eulemur_macaco            G\n"
            "Lemur_catta               G\n"
            "Nomascus_concolor         G\n"
            "Alouatta_seniculus        V\n"
            "Eulemur_rubriventer       V\n"
            "Macaca_thibetana          V\n"
            "Macaca_tonkeana           V\n"
            "Mico_humeralifer          I\n"
            "Trachypithecus_pileatus   I\n"
        )
        trait_config.write_text(
            "Aotus_trivirgatus\t1\n"
            "Ateles_geoffroyi\t1\n"
            "Cercocebus_torquatus\t1\n"
            "Eulemur_fulvus\t1\n"
            "Alouatta_seniculus\t0\n"
            "Eulemur_rubriventer\t0\n"
            "Macaca_thibetana\t0\n"
            "Macaca_tonkeana\t0\n"
        )
        resampled_traits.write_text(
            "b_1\tAotus_trivirgatus,Ateles_geoffroyi,Cercocebus_torquatus,Eulemur_fulvus\tAlouatta_seniculus,Eulemur_rubriventer,Macaca_thibetana,Macaca_tonkeana\n"
            "b_2\tAotus_trivirgatus,Ateles_geoffroyi,Cercocebus_torquatus,Eulemur_fulvus\tAlouatta_seniculus,Eulemur_rubriventer,Macaca_thibetana,Mico_humeralifer\n"
            "b_3\tAotus_trivirgatus,Ateles_geoffroyi,Cercocebus_torquatus,Eulemur_macaco\tAlouatta_seniculus,Eulemur_rubriventer,Macaca_thibetana,Macaca_tonkeana\n"
            "b_4\tEulemur_fulvus,Eulemur_macaco,Lemur_catta,Nomascus_concolor\tMacaca_thibetana,Macaca_tonkeana,Mico_humeralifer,Trachypithecus_pileatus\n"
        )
        return alignment, trait_config, resampled_traits

    def run_bootstrap(self, alignment, trait_config, resampled_traits, output, *arguments):
        command = [
            sys.executable,
            str(CT),
            "bootstrap",
            "-a", str(alignment),
            "-t", str(trait_config),
            "-s", str(resampled_traits),
            "-o", str(output),
            "--fmt", "phylip-relaxed",
        ]
        command.extend(arguments)
        return subprocess.run(
            command,
            cwd=CAASTOOLS_DIR,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

    def test_cli_writes_one_merged_species_deduplicated_event(self):
        with tempfile.TemporaryDirectory() as directory:
            alignment, trait_config, resampled_traits = self.write_compatible_inputs(directory)
            output = Path(directory) / "compatible.bootstrap.tsv"
            legacy_only_output = Path(directory) / "compatible.legacy-only.bootstrap.tsv"
            events_output = Path(directory) / "compatible.events.tsv"

            run = self.run_bootstrap(
                alignment,
                trait_config,
                resampled_traits,
                output,
                "--summarize-species", "yes",
                "--event-output", str(events_output),
            )

            self.assertEqual(run.returncode, 0, run.stdout)
            legacy_only_run = self.run_bootstrap(
                alignment,
                trait_config,
                resampled_traits,
                legacy_only_output,
            )
            self.assertEqual(legacy_only_run.returncode, 0, legacy_only_run.stdout)
            self.assertEqual(output.read_text(), legacy_only_output.read_text())
            legacy_fields = output.read_text().strip().split("\t")
            self.assertEqual(legacy_fields[1], "3")
            self.assertEqual(legacy_fields[4], "b_1,b_2,b_3")

            with events_output.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(len(rows), 1)
            event = rows[0]
            self.assertEqual(list(event), EVENT_HEADER)
            self.assertEqual(event["primary_event"], "yes")
            self.assertEqual(event["fg_amino_acids"], "A,G")
            self.assertEqual(event["fg_dominant_amino_acid"], "A=0.8")
            self.assertEqual(event["bg_amino_acids"], "I,V")
            self.assertEqual(event["bg_dominant_amino_acid"], "V=0.8")
            self.assertEqual(event["hypotheses"], "b_1,b_2,b_3")
            self.assertEqual(event["event_pattern"], "4")
            self.assertEqual(event["fg_support_count"], "5")
            self.assertEqual(event["bg_support_count"], "5")
            self.assertEqual(event["fg_observed_count"], "7")
            self.assertEqual(event["bg_observed_count"], "6")
            self.assertEqual(event["fg_discovery_count"], "7")
            self.assertEqual(event["bg_discovery_count"], "6")
            self.assertEqual(event["fg_support_amino_acid_counts"], "A=4,G=1")
            self.assertEqual(event["bg_support_amino_acid_counts"], "I=1,V=4")
            self.assertEqual(event["conflicting_species_count"], "0")
            self.assertIn("Aotus_trivirgatus=A", event["fg_support_species_amino_acids"])
            self.assertNotIn("f1=", event["fg_support_species_amino_acids"])
            self.assertIn("Fixed discovery pools: FG=7, BG=6", run.stdout)

    def test_incompatible_signatures_remain_separate_events(self):
        processed_position = SimpleNamespace(
            d={
                "f1": "A@12",
                "f2": "G@12",
                "b1": "V@12",
                "b2": "A@12",
            }
        )
        positive_hypotheses = [
            {
                "hypothesis": "b_1",
                "pattern": "1",
                "fg_amino_acids": ["A"],
                "bg_amino_acids": ["V"],
                "fg_species_amino_acids": {"f1": "A"},
                "bg_species_amino_acids": {"b1": "V"},
            },
            {
                "hypothesis": "b_2",
                "pattern": "1",
                "fg_amino_acids": ["G"],
                "bg_amino_acids": ["A"],
                "fg_species_amino_acids": {"f2": "G"},
                "bg_species_amino_acids": {"b2": "A"},
            },
        ]

        events = summarize_position_events(
            processed_position,
            positive_hypotheses,
            discovery_fg=["f1", "f2"],
            discovery_bg=["b1", "b2"],
        )

        self.assertEqual(len(events), 2)
        signatures = {
            (frozenset(event["fg_amino_acids"]), frozenset(event["bg_amino_acids"]))
            for event in events
        }
        self.assertEqual(
            signatures,
            {(frozenset({"A"}), frozenset({"V"})), (frozenset({"G"}), frozenset({"A"}))},
        )
        self.assertTrue(all(
            len(event["conflicting_fg"]) + len(event["conflicting_bg"]) == 1
            for event in events
        ))

    def test_nontransitive_compatibility_yields_all_maximal_events(self):
        processed_position = SimpleNamespace(
            d={"fA": "A@4", "fG": "G@4", "bV": "V@4", "bA": "A@4"}
        )
        positive_hypotheses = [
            {
                "hypothesis": "b_1",
                "pattern": "1",
                "fg_amino_acids": ["A"],
                "bg_amino_acids": ["V"],
                "fg_species_amino_acids": {"fA": "A"},
                "bg_species_amino_acids": {"bV": "V"},
            },
            {
                "hypothesis": "b_2",
                "pattern": "1",
                "fg_amino_acids": ["G"],
                "bg_amino_acids": ["V"],
                "fg_species_amino_acids": {"fG": "G"},
                "bg_species_amino_acids": {"bV": "V"},
            },
            {
                "hypothesis": "b_3",
                "pattern": "1",
                "fg_amino_acids": ["G"],
                "bg_amino_acids": ["A"],
                "fg_species_amino_acids": {"fG": "G"},
                "bg_species_amino_acids": {"bA": "A"},
            },
        ]

        events = summarize_position_events(
            processed_position,
            positive_hypotheses,
            discovery_fg=["fA", "fG"],
            discovery_bg=["bA", "bV"],
        )

        self.assertEqual(len(events), 2)
        event_hypotheses = {tuple(event["hypotheses"]) for event in events}
        self.assertEqual(event_hypotheses, {("b_1", "b_2"), ("b_2", "b_3")})

    def test_exact_event_pvalue_uses_fixed_group_sizes(self):
        pvalue = exact_signature_separation_pvalue(
            fg_amino_acids={"A"},
            bg_amino_acids={"V"},
            observed_fg={"f1": "A", "f2": "A"},
            observed_bg={"b1": "V", "b2": "V"},
        )
        self.assertAlmostEqual(pvalue, 1.0 / 6.0)

    def test_cross_side_species_are_rejected_only_for_event_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            alignment, trait_config, resampled_traits = self.write_compatible_inputs(directory)
            with resampled_traits.open("a") as handle:
                handle.write(
                    "b_5\tAlouatta_seniculus,Aotus_trivirgatus,Ateles_geoffroyi,Cercocebus_torquatus"
                    "\tEulemur_rubriventer,Macaca_thibetana,Macaca_tonkeana,Mico_humeralifer\n"
                )
            output = Path(directory) / "cross-side.bootstrap.tsv"

            legacy_run = self.run_bootstrap(
                alignment, trait_config, resampled_traits, output
            )
            event_run = self.run_bootstrap(
                alignment,
                trait_config,
                resampled_traits,
                output,
                "--summarize-species", "yes",
            )

            self.assertEqual(legacy_run.returncode, 0, legacy_run.stdout)
            self.assertNotEqual(event_run.returncode, 0)
            self.assertIn("fixed FG/BG membership", event_run.stdout)


if __name__ == "__main__":
    unittest.main()
