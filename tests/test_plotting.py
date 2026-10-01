from __future__ import annotations

import unittest

from hcrprobeforge.core import (
    _annotation_label_mode,
    _grouped_probe_label_modes,
    _uniform_probe_label_mode,
)


class PlotLabelPolicyTests(unittest.TestCase):
    def test_annotation_labels_are_decided_per_region(self):
        # A short 5' UTR must not force a wide CDS label into a callout lane.
        self.assertEqual(
            _annotation_label_mode(
                kind="5UTR",
                width_nt=35,
                label_width_nt=42,
                premrna_map=False,
            ),
            "outside",
        )
        self.assertEqual(
            _annotation_label_mode(
                kind="CDS",
                width_nt=1200,
                label_width_nt=30,
                premrna_map=False,
            ),
            "inside",
        )

    def test_narrow_premrna_features_are_omitted(self):
        self.assertEqual(
            _annotation_label_mode(
                kind="exon",
                width_nt=20,
                label_width_nt=45,
                premrna_map=True,
            ),
            "omit",
        )
        self.assertEqual(
            _annotation_label_mode(
                kind="intron",
                width_nt=900,
                label_width_nt=45,
                premrna_map=True,
            ),
            "inside",
        )

    def test_similar_probe_boxes_use_one_inside_mode(self):
        self.assertEqual(
            _uniform_probe_label_mode([50.0, 52.0, 48.0], [7.0, 7.0, 7.0]),
            "inside",
        )

    def test_similar_probe_boxes_use_one_outside_mode_when_none_fit(self):
        self.assertEqual(
            _uniform_probe_label_mode([10.0, 10.5, 11.0], [7.0, 7.0, 7.0]),
            "outside",
        )

    def test_materially_different_probe_boxes_allow_individual_fallback(self):
        self.assertIsNone(
            _uniform_probe_label_mode([10.0, 20.0], [7.0, 7.0])
        )

    def test_grouped_probe_modes_keep_equal_sized_boxes_consistent(self):
        modes = _grouped_probe_label_modes(
            [52.0, 52.0, 52.0, 36.0],
            [7.0, 7.0, 7.0, 7.0],
        )
        self.assertEqual(len(set(modes[:3])), 1)

    def test_grouped_probe_modes_can_differ_for_different_sizes(self):
        modes = _grouped_probe_label_modes(
            [52.0, 52.0, 36.0],
            [7.0, 7.0, 40.0],
        )
        self.assertEqual(modes, ["inside", "inside", "outside"])

    def test_adjacent_short_tiers_do_not_split_similar_50_and_52_nt_boxes(self):
        # This reproduces the previous phox2bb map: 38/40-nt boxes were
        # grouped with 50-nt boxes, while 52-nt boxes were separated. The
        # 50- and 52-nt boxes must therefore share one placement mode.
        modes = _grouped_probe_label_modes(
            [13.73, 14.45, 18.07, 18.79],
            [9.0, 7.0, 7.0, 7.0],
        )
        self.assertEqual(modes, ["outside", "outside", "inside", "inside"])


if __name__ == "__main__":
    unittest.main()
