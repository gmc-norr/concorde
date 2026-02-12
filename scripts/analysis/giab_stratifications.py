"""GIAB stratification BED file loader.

Discovers and maps GIAB benchmark stratification BED files from a
standard GIAB directory layout to Concorde's region_beds configuration.

GIAB stratification directories typically follow this structure:
    <bed_dir>/
        LowComplexity/
            GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz
            ...
        SegmentalDuplications/
            GRCh38_segdups.bed.gz
            ...
        GCcontent/
            GRCh38_gc15_slop50.bed.gz
            ...
        Mappability/
            GRCh38_lowmappabilityall.bed.gz
            ...

This module maps these to our internal dimension names:
    LowComplexity -> low_complexity
    SegmentalDuplications -> segmental_duplications (segdup)
    GCcontent -> gc_content
    Mappability -> mappability
"""

import logging
from pathlib import Path

log = logging.getLogger(__name__)

# Map GIAB directory names to our internal dimension keys
_GIAB_CATEGORY_MAP = {
    "LowComplexity": "low_complexity",
    "SegmentalDuplications": "segmental_duplications",
    "GCcontent": "gc_content",
    "Mappability": "mappability",
}


class GIABStratificationLoader:
    """Load GIAB benchmark stratification BED files.

    Args:
        bed_dir: Root directory of GIAB stratification BEDs
        categories: List of GIAB category names to load (None = all available)
    """

    def __init__(self, bed_dir, categories=None):
        self.bed_dir = Path(bed_dir) if bed_dir else None
        self.categories = categories

    def discover_bed_files(self):
        """Discover BED files organized by GIAB category.

        Returns:
            Dict mapping category name to list of BED file paths.
        """
        if not self.bed_dir or not self.bed_dir.is_dir():
            log.warning("GIAB bed_dir not found or not configured: %s", self.bed_dir)
            return {}

        discovered = {}
        for category_dir in sorted(self.bed_dir.iterdir()):
            if not category_dir.is_dir():
                continue

            category_name = category_dir.name
            if self.categories and category_name not in self.categories:
                continue

            bed_files = sorted(category_dir.glob("*.bed.gz"))
            if bed_files:
                discovered[category_name] = [str(f) for f in bed_files]
                log.info(
                    "Discovered %d BED files in %s",
                    len(bed_files),
                    category_name,
                )

        return discovered

    def get_region_beds_config(self):
        """Generate region_beds config dict from GIAB BEDs.

        Uses the first BED file in each recognized category as the default.
        For categories with multiple BED files, the first (alphabetically)
        is selected.

        Returns:
            Dict suitable for merging into config region_beds section.
        """
        discovered = self.discover_bed_files()
        config = {}

        for giab_name, internal_name in _GIAB_CATEGORY_MAP.items():
            beds = discovered.get(giab_name, [])
            if beds:
                config[internal_name] = beds[0]
                log.info("Mapped GIAB %s -> %s: %s", giab_name, internal_name, beds[0])

        return config

    def get_all_bed_files(self):
        """Get all discovered BED files as a flat dict.

        Returns:
            Dict mapping (category, filename) to path strings.
        """
        discovered = self.discover_bed_files()
        result = {}
        for category, paths in discovered.items():
            for p in paths:
                filename = Path(p).name
                result[(category, filename)] = p
        return result
